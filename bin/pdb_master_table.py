import sys
import glob
from os import path
from argparse import ArgumentParser
import re
import warnings

import pandas

from Bio import SeqIO
from Bio import PDB
from Bio.Seq import Seq
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
# argparser.add_argument('--alnroot', metavar='aln_root', type=str, required=True)

args = argparser.parse_args()


pdb_dir = args.pdbdir
slr_root = args.slrroot

# Ensembl protein sequences
# ens_seqs = dict()
# with open("data/Homo_sapiens.GRCh37.58.pep.all.fa", "rU") as ensembl:
#     ens_fasta = SeqIO.parse(ensembl, "fasta")
#     for record in ens_fasta:
#         ens_seqs[record.id] = record.seq

# uniprot_sprot.fasta comes from:
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/
print >>sys.stderr, "Loading Uniprot sequences...",
seqs = dict()
with open("data/uniprot_sprot.fasta", "rU") as uniprot:
    fasta = SeqIO.parse(uniprot, "fasta")
    for record in fasta:
        up = re_uniprotid.match(record.id).groups()[0]
        seqs[up] = record.seq
print >>sys.stderr, "done."

pdb_map_file = open(args.pdbmap)
outfile = open(args.outfile, 'w')

p = PDB.PDBParser()
# p = PDB.MMCIFParser()

for l in pdb_map_file:
# for l in [ "ENSP00000360777	2xsw	B	Q9NRR6		2	350	275	623	275	623\n" ]:
    f = l.rstrip().split('\t')
    ens, pdb_name, chain_name, uniprot = f[:4]
    print uniprot, pdb_name, chain_name

    # pdb_fn = path.join(pdb_dir, 'pdb' + pdb_name + '.cif')
    pdb_fn = path.join(pdb_dir, 'pdb' + pdb_name + '.ent')
    pdb = p.get_structure(ens, pdb_fn)

    # TODO Access the header SEQRES information?
    pdb_begin_id, pdb_begin_ins = parse_coord(f[6])
    pdb_end_id, pdb_end_ins = parse_coord(f[7])
    uniprot_begin, uniprot_end = int(f[8])-1, int(f[9])

    # TODO: Make sure we only have X-Ray structures, not NMR
    print "# of models", len(pdb)
    model = pdb[0]
    dssp = DSSP(model, pdb_fn) 
    chain = pdb[0][chain_name]

    sites_fn_glob = glob.glob(path.join(slr_root, args.clade, '*', ens+'*'))
    if not len(sites_fn_glob):
        print "Missing results file for", ens
        continue

    sites_fn = sites_fn_glob[0]

    pdb_chain = list(chain)
    print pdb_begin_id, pdb_begin_ins
    # print pdb_begin_id, pdb_chain[0].id, pdb_end_id, pdb_chain[-1].id
    # TODO Sergio mentioned that in cases where residue's orientation couldn't be determined exactly,
    # there may be multiple entries with a the same number and different letters
    # -> I could explicitly check for that in the code?

    # TODO Make a function out of this whole mapping business?

    # TODO Instead of looping we could try to remove all non-AA residues from the chain?
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

    print pdb_begin_id, found_begin_id, pdb_end_id, found_end_id

    # The indexing offset between Ensembl and PDB
    # TODO Do we need an extra 1 to account for 1- vs. 0-based indexing?
    offset = int(f[9]) - int(f[7]) - 1
    # print "Offset", offset

    # site_map = {}
    site_map = []
    # sites = pandas.read_fwf(open(sites_fn), colspecs=colspecs, header=False, comment="\n")
    # for l in sites.iterrows():
    #   print l[0], l[1]
    #   site_map[l[0]+offset] = l[1]
    sites = open(sites_fn)
    sites.readline()
    for i, l in enumerate(sites):
        f = l.rstrip().split('\t')
        # I removed offset when rewriting the matching process to use pairwise alignment
        # site_map[int(f[0])+offset] = float(f[3])
        # TODO Does this fix the mapping issue?
        # site_map[i] = float(f[3])
        site_map.append(float(f[3]))
        # site_map[int(f[0])] = float(f[3])

    # print "Len of the uniprot seq", len(seqs[uniprot])
    # print ''.join([ SeqUtils.seq1(r.resname) for r in pdb_chain[found_begin_i:found_end_i+1] ])
    # print found_begin_i+offset, found_end_i+offset
    # print seqs[uniprot][uniprot_begin:uniprot_end]

    # for s in SeqIO.parse(pdb_fn, 'pdb-seqres'):
    #    if s.id == pdb_name.upper() + ':' + chain_name:
    #        pdb_seq = s.seq

    pdb_seq = ''.join([ SeqUtils.seq1(r.resname) for r in pdb_chain[found_begin_i:found_end_i+1] ])

    alignment = align.globalxs(seqs[uniprot][uniprot_begin:uniprot_end], pdb_seq, -5, -1)

    if not len(alignment):
        continue
    # TODO I cut corners here by only taking into acc. the first alignment
    # print format_alignment(*alignment[0])
    n_matched = 0
    
    pdb_vals = [None] * len(pdb_seq)
    ref = alignment[0][0]
    pdb_target = alignment[0][1]
    ref_i, pdb_i = 0, 0
    for i, r in enumerate(ref):
        if r != '-':
            if pdb_target[i] != '-':
                # pdb_vals[pdb_i] = site_map.get(ref_i+uniprot_begin, float('nan'))
                pdb_vals[pdb_i] = site_map[ref_i+uniprot_begin]
                n_matched += 1
            ref_i += 1

        if pdb_target[i] != '-':
            pdb_i += 1

    # print pdb_vals
    n_mismatched = uniprot_end - uniprot_begin - n_matched
    print n_mismatched, "mismatched out of", found_end_i - found_begin_i

    for i, r in enumerate(pdb_chain[found_begin_i:found_end_i+1]):
        # omega = site_map.get(i, 'nan')
        omega = pdb_vals[i]
        
        i = i + found_begin_id
        # print SeqUtils.seq1(r.resname), seqs[uniprot][i+offset]
        # TODO Now make sure the remaining bit is all aas?
        # Also taking into account modified aas?
        if not (PDB.is_aa(r) or r.resname in ['OLD']):
            print "Dodgy PDB residue", r
            continue
            # raise ValueError("Not an amino-acid")

        # TODO Compare to the sequence from Uniprot?
        # print SeqUtils.seq1(r.resname), i+offset, seqs[uniprot][i+offset]
        # if SeqUtils.seq1(r.resname) != seqs[uniprot][i+offset]:
        #     n_mismatched += 1
            # raise ValueError("Residue doesn't agree with reference")

        # TODO Add a check for residue number agreement
        # and/or check sequence context from the ptmfun file
        try:
            ann = dssp[(chain_name, r.id)]
        except KeyError:
            continue

        print >>outfile, '\t'.join((str(i) for i in (ens, pdb_name, r.id[1],  ann[1],
                                                     ann[3], omega)))
