import re
import sys
import glob
from os import path
from argparse import ArgumentParser
import utils

from Bio import SeqIO
from Bio import AlignIO
from Bio import SeqUtils
from Bio import PDB

import numpy as np
from sklearn.cluster import AffinityPropagation

import brewer2mpl

re_resid = re.compile("(-?[0-9]+)([A-Z]*)")

def parse_coord(coord):
    n, ic = re_resid.match(coord).groups()
    if ic == '':
        ic = ' '

    return int(n), ic

p = PDB.PDBParser()

argparser = ArgumentParser()
argparser.add_argument("--pdbmap", metavar="pdb_map", type=str, required=True)
argparser.add_argument("--ensid", metavar="ens_id", type=str, required=True)
argparser.add_argument("--pdbfile", metavar="pdb_file", type=str, required=True)
argparser.add_argument("--alnfile", metavar="aln_file", type=str, required=True)
argparser.add_argument("--mapfile", metavar="map_file", type=str, required=True)
argparser.add_argument("--colorfile", metavar="color_file", type=str, required=True)

args = argparser.parse_args()

for l in open(args.pdbmap):
    f = l.rstrip().split()
    if f[0] == args.ensid:
        chain_id = f[2]
        pdb_begin_id, pdb_begin_ins = parse_coord(f[6])
        pdb_end_id, pdb_end_ins = parse_coord(f[7])
        break
else:
    raise ValueError("Ensembl ID not found")

omega_map = []
# Omega   lower   upper   LRT_Stat        Pval    Adj.Pval
for l in open(args.mapfile):
    f = l.rstrip().split()
    omega_map.append(f)

seqs = list(AlignIO.read(args.alnfile, 'fasta'))

pdb_seq = seqs[0]
ens_seq = seqs[1]

pdb_id = pdb_seq.name
ens_id = ens_seq.name

pdb = p.get_structure(pdb_id, args.pdbfile)
chain = list(pdb[0][chain_id])

found_begin_id, found_begin_i, found_begin_id, found_end_i = \
utils.parse_chain(chain, pdb_begin_id, pdb_begin_ins, pdb_end_id, pdb_end_ins)

chain = chain[found_begin_i:found_end_i+1]
chain_seq = ''.join([ SeqUtils.seq1(r.resname) for r in chain if SeqUtils.seq1(r.resname) != 'X' ])
# assert chain_seq == str(pdb_seq.seq.ungap('-'))
print chain_seq
print str(pdb_seq.seq.ungap('-'))

r_coords = []
residues = []
for i, r in enumerate(chain):
    if float(omega_map[i][0]) > 1:
        r_coords.append(r['CA'].get_coord())
        residues.append(''.join((str(_id) for _id in r.id)))

r_coords = np.array(r_coords)
af = AffinityPropagation().fit(r_coords)
labels = af.labels_

palette = brewer2mpl.get_map('Spectral', 'Diverging', max(labels)+1)
colorfile = open(args.colorfile, 'w')

for i, r in enumerate(residues):
    l = labels[i]
    print >>colorfile, '\t'.join(str(it) for it in [r] + palette.colors[l])
