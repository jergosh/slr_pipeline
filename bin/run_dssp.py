import sys
from os import path
import glob
import re
from argparse import ArgumentParser

from Bio import SeqIO, SeqUtils
from Bio import PDB
from Bio.PDB.DSSP import DSSP
from Bio.PDB.PDBExceptions import PDBException
from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62

re_uniprotid = re.compile(".*\|(.*)\|.*")
re_resid = re.compile("(-?[0-9]+)([A-Z]*)")

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

dssp_missing = 0
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
    with open("data/uniprot_sprot.fasta", "rU") as uniprot:
        fasta = SeqIO.parse(uniprot, "fasta")
        for record in fasta:
            up = re_uniprotid.match(record.id).groups()[0]
            seqs[up] = record.seq
    print >>sys.stderr, "done."


    previous_pdb_id = None
    dssp_result = None

    sifts_map_file = open(args.siftsmap)
    outfile = open(args.outfile, 'w')

    for l in sifts_map_file:
        ens, uniprot, uniprot_coord, pdb_id, pdb_chain, res_id = l.rstrip().split('\t')
        # We can make the decision about a 
        if pdb_id == previous_pdb_id: # We're in the middle of a gene
            if dssp_result == None:
                # TODO do we need to set some state here?
                continue            

        else: # New gene
            print "Processing", ens
            slr_fn_glob = glob.glob(path.join(args.slrroot, args.clade, '*', ens+'*'))
            if not len(slr_fn_glob):
                print "Missing results file for", ens
                continue

            slr_fn = slr_fn_glob[0]
            sites = open(slr_fn)
            sites.readline()
            site_map = []
            for i, l in enumerate(sites):
                f = l.rstrip().split('\t')
                site_map.append(float(f[3]))

            pdb_fn = path.join(args.pdbdir, 'pdb'+pdb_id+'.ent')
            try:
                pdb = p.get_structure(pdb_id, pdb_fn)
            except IOError, e:
                print >>sys.stderr, "Missing PDB file for", pdb_id
                dssp_result = None
                continue
            # TODO Once again, not entirely clear what to do about NMR structures 
            model = pdb[0]
            previous_pdb_id = pdb_id

            alignment = None
            if str(ens_seqs[ens]) != str(seqs[uniprot]):
                alignment = align.localxs(seqs[uniprot], ens_seqs[ens], -10, -0.5)[0]
                pos_map = aln_to_map(alignment)
            else:
                pos_map = None

            try:
                dssp_result = DSSP(model, pdb_fn)
            except PDBException, e:
                dssp_result = None
                continue
        
        uniprot_coord = int(uniprot_coord)-1
        if alignment is None:
            coord = uniprot_coord
        else:
            try:
                coord = pos_map[uniprot_coord]
            except KeyError, e:
                print "Unaligned residue", uniprot_coord, "in", uniprot
                continue

        omega = site_map[coord]
        # print omega

        try:
            dssp_info = dssp_result[(pdb_chain, parse_coord(res_id))]
            print >>outfile, '\t'.join([ ens, str(coord), uniprot, str(uniprot_coord), pdb_id, pdb_chain, res_id,
                                         str(dssp_info[1]), str(dssp_info[3]), str(omega) ])
        except KeyError, e:
            dssp_missing += 1

            
    print >>sys.stderr, "Missing from DSSP", dssp_missing

