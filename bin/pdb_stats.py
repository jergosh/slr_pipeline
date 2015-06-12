# FIXME No longer gives the right statistic for overall coverage as the map file
# can contain multiple entries for each protein identifier.
import sys
from argparse import ArgumentParser
import urllib, urllib2
import glob
from os import path

from Bio import SeqIO

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

argparser = ArgumentParser()

# argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--pep', metavar='pep_file', type=str, required=True)
argparser.add_argument('--mapfile', metavar='sifts_file', type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

args = argparser.parse_args()

# How to extract the CDS lengths?
# the .pep file?
seq_dict = SeqIO.to_dict(SeqIO.parse(args.pep, 'fasta'))

ens_sum, pdb_sum = 0, 0
coverages, coverages_rel = [], []
for l in open(args.mapfile):
    f = l.rstrip().split('\t')
    coverage = float(int(f[9])-int(f[8]))

    coverages.append(coverage)
    coverages_rel.append(min(1, coverage/len(seq_dict[f[0]].seq)))

    pdb_sum += coverage
    ens_sum += len(seq_dict[f[0]].seq)

pages = PdfPages(args.outfile)

fig = plt.figure()
n, bins, patches = plt.hist(coverages, 50, normed=0, histtype='stepfilled')
plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
plt.title("PDB coverage")
plt.xlabel("# residues")
plt.ylabel("Count")
pages.savefig(fig, bbox_inches=0)

fig = plt.figure()
n, bins, patches = plt.hist(coverages_rel, 50, normed=0, histtype='stepfilled')
plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
plt.title("PDB coverage")
plt.xlabel("Fraction")
plt.ylabel("Count")
pages.savefig(fig, bbox_inches=0)
pages.close()

print float(pdb_sum)/ens_sum
