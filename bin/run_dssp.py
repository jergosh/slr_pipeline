import sys
from os import path
import glob
import re
from argparse import ArgumentParser

import pandas
from Bio import SeqIO, SeqUtils
from Bio import PDB
from Bio.PDB.DSSP import DSSP
from Bio.PDB.PDBExceptions import PDBException
from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62

import slr

re_uniprotid = re.compile(".*\|(.*)\|.*")
re_resid = re.compile("(-?[0-9]+)([A-Z]*)")

pdbl=PDB.PDBList()

def parse_coord(coord):
    n, ic = re_resid.match(coord).groups()
    if ic == '':
        ic = ' '

    return ' ', int(n), ic

def aln_to_map(aln):
    "Assuming reference is the first one in the alignment"
    pos_map = {}
    
    seqA, seqB, score, begin, end = aln

    iA, iB = 0, 0
    for i, p in enumerate(seqA):
        if p != '-':
            if seqB[i] != '-':
                if i >= begin and i < end:
                    pos_map[iA] = iB

            iA += 1

        if seqB[i] != '-':
            iB += 1

    return pos_map

p = PDB.PDBParser(QUIET=True)

argparser = ArgumentParser()

argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--siftsmap', metavar='sifts_map', type=str, required=True)
argparser.add_argument('--refgenome', metavar='ref_genome', type=str, required=True)
argparser.add_argument('--pdbdir', metavar='pdb_dir', type=str, required=True)
argparser.add_argument('--slrroot', metavar='slr_root', type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

def process_slr(df, slrroot, clade, pdbdir):
    dssp_matched, dssp_missing = 0, 0

    ens = df['stable_id'].iloc[0]
    uniprot = df['uniprot_id'].iloc[0]
    pdb_id = df['pdb_id'].iloc[0]
    assert len(set(df['pdb_chain'])) == 1
    pdb_chain = df['pdb_chain'].iloc[0]
    print "Processing", ens, pdb_id, pdb_chain

    slr_fn_glob = glob.glob(path.join(slrroot, clade, '*', ens+'*'))
    if not len(slr_fn_glob):
        print "Missing results file for", ens
        return

    slr_fn = slr_fn_glob[0]
    sites = open(slr_fn)
    sites.readline()
    site_map = []
    for i, l in enumerate(sites):
        f = l.rstrip().split('\t')
        site_map.append(float(f[3]))

    if len(site_map) != len(ens_seqs[ens]):
        assert len(site_map) == len(ens_seqs[ens])+1

    # pdb_fn = path.join(pdbdir, 'pdb'+pdb_id+'.ent')
    pdb_fn = path.join(pdbdir, pdb_id)
    try:
        if not path.exists(pdb_fn):
            print "Fetching", pdb_id
            pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=args.pdbdir)
            print pdb_file

        pdb = p.get_structure(pdb_id, pdb_fn)
    except IOError, e:
        print >>sys.stderr, "Missing PDB file for", pdb_id
        dssp_result = None
        return

    # TODO Once again, not entirely clear what to do about NMR structures
    if len(pdb) > 1:
        print >>sys.stderr, "Multiple models", len(pdb)
        return

    if not seqs.get(uniprot):
        print >>sys.stderr, "Uniprot sequence missing"
        return

    model = pdb[0]
    alignment = None
    if str(ens_seqs[ens]) != str(seqs[uniprot]):
        alignment = align.localxs(seqs[uniprot], ens_seqs[ens], -10, -0.5)[0]
        pos_map = aln_to_map(alignment)
    else:
        pos_map = None

    try:
        dssp_result = DSSP(model, pdb_fn)
    except PDBException, e:
        print >>sys.stderr, "PDB file missing?"
        # TODO Consider going through with the mapping
        return
    except Exception, e:
        print >>sys.stderr, "DSSP result missing"
        # TODO Consider going through with the mapping
        return

    for index, row in df.iterrows():
        uniprot_coord = row['uniprot_coord']-1
        if alignment is None:
            coord = uniprot_coord
        else:
            try:
                coord = pos_map[uniprot_coord]
            except KeyError, e:
                print "Unaligned residue", uniprot_coord, "in", uniprot
                continue

        omega = site_map[coord]

        try:
            res_id = str(row['res_id'])
            dssp_info = dssp_result[(pdb_chain, parse_coord(res_id))]
            print >>outfile, '\t'.join([ ens,
                                         str(coord),
                                         uniprot,
                                         str(uniprot_coord),
                                         pdb_id, pdb_chain, res_id,
                                         str(dssp_info[2]),
                                         str(dssp_info[3]),
                                         str(omega) ])
        except KeyError, e:
            print "DSSP KeyError", res_id
            dssp_missing += 1
        else:
            dssp_matched += 1

    print dssp_matched, dssp_missing


if __name__ == "__main__":
    args = argparser.parse_args()

    print >>sys.stderr, "Loading reference genome...",
    ens_seqs = dict()
    with open(args.refgenome, "r") as uniprot:
        fasta = SeqIO.parse(uniprot, "fasta")
        for record in fasta:
            ens_seqs[record.name] = record.seq
    print >>sys.stderr, "done."

    print >>sys.stderr, "Loading Uniprot sequences...",
    seqs = dict()
    with open("/nfs/research1/goldman/gregs/slr_pipeline/data/uniprot_sprot.fasta", "rU") as uniprot:
        fasta = SeqIO.parse(uniprot, "fasta")
        for record in fasta:
            up = re_uniprotid.match(record.id).groups()[0]
            seqs[up] = record.seq
    print >>sys.stderr, "done."

    outfile = open(args.outfile, 'w')

    sifts_map = pandas.read_csv(open(args.siftsmap), header=None, comment="\n", sep="\t")
    sifts_map.columns = [ "stable_id", "uniprot_id", "uniprot_coord", "pdb_id", "pdb_chain", "res_id" ]
    sifts_map.groupby([ "stable_id", "pdb_id", "pdb_chain" ]).apply(process_slr, args.slrroot, args.clade, args.pdbdir)

    # print >>sys.stderr, "Missing from DSSP", dssp_missing
