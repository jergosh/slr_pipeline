import re
import sys
import glob
from os import path
from argparse import ArgumentParser
import random
import utils
from collections import defaultdict
from pprint import pprint

from Bio import SeqIO
from Bio import AlignIO
from Bio import SeqUtils
from Bio import PDB

import numpy as np
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import KMeans

palette = [[99,204,82],
[189,89,212],
[211,77,41],
[59,64,57],
[202,207,196],
[213,168,59],
[199,150,200],
[85,110,42],
[83,59,105],
[189,208,130],
[106,218,182],
[196,70,141],
[104,116,205],
[185,133,128],
[193,218,60],
[90,138,114],
[102,42,42],
[107,164,195],
[207,74,91],
[174,116,62]]

re_resid = re.compile("(-?[0-9]+)([A-Z]*)")

def isnan(f):
    return not (f == f)

def parse_coord(coord):
    n, ic = re_resid.match(coord).groups()
    if ic == '':
        ic = ' '

    return int(n), ic

d = lambda: defaultdict(int)
cl_found = defaultdict(d)
def collapse_clusters(labels):
    clusters = defaultdict(int)
    cluster_sizes = defaultdict(int)
    for l in labels:
        if isnan(l):
            print >>sys.stderr, labels
            return {}
        clusters[l] += 1

    for cl_size in clusters.values():
        cluster_sizes[cl_size] += 1

    return cluster_sizes

def process_clusters(labels):
    cluster_sizes = collapse_clusters(labels)

    for cl_size in cluster_sizes:
        cl_found[cl_size][cluster_sizes[cl_size]] += 1

def get_pvals(cluster_found, niter):
    for n in sorted(cl_found, reverse=True):
        for k in sorted(cl_found[n], reverse=True):
            if n > 1:
                cl_found[n-1][k] += cl_found[n][k]
            if k > 1:
                cl_found[n][k-1] += cl_found[n][k]

            cl_found[n][k] = min(cl_found[n][k]/float(niter), 1.0)

p = PDB.PDBParser(QUIET=True)

argparser = ArgumentParser()
argparser.add_argument("--pdbmap", metavar="pdb_map", type=str, required=True)
argparser.add_argument("--ensid", metavar="ens_id", type=str, required=True)
argparser.add_argument("--pdbfile", metavar="pdb_file", type=str, required=True)
argparser.add_argument("--alnfile", metavar="aln_file", type=str, required=True)
argparser.add_argument("--mapfile", metavar="map_file", type=str, required=True)
argparser.add_argument("--colorfile", metavar="color_file", type=str, required=True)

argparser.add_argument("--niter", metavar="n_iter", type=int, required=False, default=0)
argparser.add_argument("--preference", metavar="preference", type=int, required=False, default=-50)


args = argparser.parse_args()

for l in open(args.pdbmap):
    f = l.rstrip().split()
    if f[0] == args.ensid:
        chain_id = f[2]
        pdb_begin_id, pdb_begin_ins = parse_coord(f[6])
        pdb_end_id, pdb_end_ins = parse_coord(f[7])
        break
else:
    raise ValueError("Ensembl ID " + args.ensid + " not found")

omega_map = []
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
chain_seq = ''.join([ SeqUtils.seq1(r.resname) for r in chain ])
# chain_seq = ''.join([ SeqUtils.seq1(r.resname) for r in chain if SeqUtils.seq1(r.resname) != 'X' ])
assert chain_seq == str(pdb_seq.seq.ungap('-'))

r_coords = []
residues = []
# Omega   lower   upper   LRT_Stat        Pval    Adj.Pval
for i, r in enumerate(chain):
    if float(omega_map[i][0]) > 1:
        r_coords.append(r['CA'].get_coord())
        residues.append(''.join((str(_id) for _id in r.id)))

if len(r_coords) < 2:
    print >>sys.stderr, "<= 1 residue with omega > 1. Exiting"
    sys.exit(0)

r_coords = np.array(r_coords)
af = AffinityPropagation(preference=args.preference, max_iter=500).fit(r_coords)
# af = KMeans(n_clusters=5).fit(r_coords)
labels = af.labels_

# palette = brewer2mpl.get_map('Spectral', 'Diverging', max(labels)+1)
if max(labels) < 20:
    colorfile = open(args.colorfile, 'w')

    for i, r in enumerate(residues):
        l = labels[i]
        print >>colorfile, '\t'.join(str(it) for it in [r] + palette[l])

print >>sys.stderr, "Running simulations...",
i = 0
while i < args.niter:
    sim_residues = random.sample(chain, len(residues)) 
    sim_coords = np.array([ r['CA'].get_coord() for r in sim_residues ])

    af = AffinityPropagation(preference=args.preference, max_iter=500).fit(sim_coords)
    sim_labels = af.labels_

    if any([ isnan(it) for it in sim_labels ]):
        # print >>sys.stderr, "NaNs in sim_labels"
        continue

    process_clusters(sim_labels)
    # af = KMeans(n_clusters=5).fit(sim_coords)
    i += 1

    
    colorfile = open(args.colorfile+'_'+str(i), 'w')

    #for r_i, r in enumerate(sim_residues):
    #    l = sim_labels[r_i]
    #    print >>colorfile, '\t'.join(str(it) for it in [''.join( [ str(it) for it in r.id] ) ] + palette[l])

print >>sys.stderr, "done."
get_pvals(cl_found, args.niter)

# for d in cl_found:
#     print d
#     print "\t", cl_found[d]

clusters = collapse_clusters(labels)
combined_pval = 1.0
for n in sorted(clusters):
    k = clusters[n]
    combined_pval *= cl_found[n][k]
    # print '{}\t{}\t{}\t{}'.format(args.ensid, n, k, cl_found[n][k])
print combined_pval

