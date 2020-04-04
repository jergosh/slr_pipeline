import sys
import glob
from os import path
from argparse import ArgumentParser
import re
import warnings
import utils

import pandas

from Bio import SeqIO
from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils
from Bio.PDB.DSSP import DSSP

from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment
# Suppress warnings from PDBParser
warnings.simplefilter("ignore")

from slr import *

def parse_coord(coord):
    n, ic = re_resid.match(coord).groups()
    if ic == '':
        ic = ' '

    return int(n), ic

re_resid = re.compile("(-?[0-9]+)([A-Z]*)")
re_uniprotid = re.compile(".*\|(.*)\|.*")

argparser = ArgumentParser()

argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--pdbmap', metavar='pdb_map', type=str, required=True)
argparser.add_argument('--pdbdir', metavar='pdb_dir', type=str, required=True)
argparser.add_argument('--slrroot', metavar='slr_root', type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)
argparser.add_argument('--mapdir', metavar='map_dir', type=str, required=True)
argparser.add_argument('--seqref', metavar='ref_fasta', type=str,
                       default="data/ens/73/pep/Homo_sapiens.GRCh37.73.pep.all.fa")
# argparser.add_argument('--alnroot', metavar='aln_root', type=str, required=True)

args = argparser.parse_args()

pdb_dir = args.pdbdir
slr_root = args.slrroot
map_dir = args.mapdir

utils.check_dir(map_dir)

print >>sys.stderr, "Loading Ensembl sequences...",
seqs = dict()
with open(args.seqref, "rU") as ensembl:
    fasta = SeqIO.parse(ensembl, "fasta")
    for record in fasta:
        ens = record.id
        seqs[ens] = record.seq
print >>sys.stderr, "done."

best_pdb_file = open(args.pdbmap)
outfile = open(args.outfile, 'w')

p = PDB.PDBParser()


for l in best_pdb_file:
# for l in [ "ENSP00000360777	2xsw	B	Q9NRR6		2	350	275	623	275	623\n" ]:
    f = l.rstrip().split('\t')
    ens, pdb_name, chain_name, uniprot = f[:4]
    print uniprot, pdb_name, chain_name

    pdb_fn = path.join(pdb_dir, 'pdb' + pdb_name + '.ent')
    pdb = p.get_structure(ens, pdb_fn)

    # TODO Access the header SEQRES information?
    pdb_begin_id, pdb_begin_ins = parse_coord(f[6])
    pdb_end_id, pdb_end_ins = parse_coord(f[7])

    # TODO: Make sure we only have X-Ray structures, not NMR
    model = pdb[0]
    dssp = DSSP(model, pdb_fn) 
    chain = pdb[0][chain_name]

    sites_fn_glob = glob.glob(path.join(slr_root, args.clade, '*', ens+'*'))
    if not len(sites_fn_glob):
        print "Missing results file for", ens
        continue

    sites_fn = sites_fn_glob[0]

    pdb_chain = [ r for r in chain]

    for i, r in enumerate(pdb_chain):
        if r.id[1] == pdb_begin_id and r.id[2] == pdb_begin_ins:
            found_begin_id = r.id[1]
            found_begin_i = i
            break
        elif r.id[1] > pdb_begin_id:
            found_begin_id = r.id[1]
            found_begin_i = i
            break
    else:
        print "Residue not found"
        continue
        # raise IndexError("Residue not found")

    for i, r in reversed(list(enumerate(pdb_chain))):
        if r.id[1:] == (pdb_end_id, pdb_end_ins):
            found_end_id = r.id[1]
            found_end_i = i
            break
        # We need to check if the residue is an aa because some chains contain misnumbered HETATMs 
        elif r.id[1] < pdb_end_id and r.id[0] == ' ' and PDB.is_aa(r):
            found_end_id = r.id[1]
            found_end_i = i
            break
    else:
        print "Residue not found"
        continue
        # raise IndexError("Residue not found")

    site_map = []

    sites = open(sites_fn)
    sites.readline()
    for i, l in enumerate(sites):
        f = l.rstrip().split('\t')
        site_map.append(f[3:9])

    ref_seq = seqs[ens]
    pdb_seq = ''.join([ SeqUtils.seq1(r.resname) for r in pdb_chain[found_begin_i:found_end_i+1] ])

    alignment = align.globalxs(ref_seq, pdb_seq, -5, -1)
    if not len(alignment):
        continue

    # print format_alignment(*alignment[0])

    n_matched = 0
    
    pdb_vals = [ ['NA'] * 6 ] * len(pdb_seq)
    ref = alignment[0][0]
    pdb_target = alignment[0][1]
    ref_i, pdb_i = 0, 0
    for i, r in enumerate(ref):
        if r != '-':
            if pdb_target[i] != '-':
                pdb_vals[pdb_i] = site_map[ref_i]
                n_matched += 1
            ref_i += 1

        if pdb_target[i] != '-':
            pdb_i += 1

    # print pdb_vals
    n_mismatched = len(pdb_seq) - n_matched
    print n_mismatched, "mismatched out of", found_end_i - found_begin_i

    # Output what we've figured out
    outdir = path.join(map_dir, ens[-2:])
    utils.check_dir(outdir)
    outfile_fa = open(path.join(outdir, ens+'.fa'), 'w')
    fa_seqs = [ SeqRecord(id=pdb_name, seq=Seq(pdb_target), description=""),
             SeqRecord(id=ens, seq=ref, description="") ]
    SeqIO.write(fa_seqs, outfile_fa, 'fasta')

    outfile_vals = open(path.join(outdir, ens+'.tab'), 'w')
    for v in pdb_vals:
        print >>outfile_vals, '\t'.join(v)

    for i, r in enumerate(pdb_chain[found_begin_i:found_end_i+1]):
        v = pdb_vals[i]
        
        if not (PDB.is_aa(r) or r.resname in ['OLD']):
            print "Dodgy PDB residue", r
            continue

        try:
            ann = dssp[(chain_name, r.id)]
        except KeyError:
            continue

        print >>outfile, '\t'.join((str(i) for i in (ens, pdb_name, r.id[1],  ann[1],
                                                     ann[3], '\t'.join(v))))
