from glob import glob
import os
from os import path
import itertools
import re
import sys
import copy
import argparse

import pandas
from Bio import AlignIO
from Bio import Alphabet
from Bio.Alphabet import IUPAC

import numpy as np

states = Alphabet.Gapped(IUPAC.IUPACProtein).letters

# Store counts in an array for speed?
# Presumably faster than a defaultdict
def counts_from_col(col):
    counts = np.zeros(21)
    for i in col:
        counts[states.find(i)] += 1

    assert counts.sum() == len(col)

    return counts

def entropy(counts):
    counts = counts[:20]
    total = counts.sum()
    e = 0

    for s in counts:
        freq = 1.0*s/total
        if freq != 0:
            e -= freq * np.log(freq)

    return e

def column_stats(counts):
    return counts[20], 1.0*counts[20]/counts.sum()

argparser = argparse.ArgumentParser()

argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--slrroot', metavar='slr_root', type=str, required=True)
argparser.add_argument('--alnroot', metavar='aln_root', type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

args = argparser.parse_args()

outfile = open(args.outfile, 'w')

for aln_fn in glob(path.join(args.alnroot, args.clade, "*", "*_prank.best.fas")):
    basename = path.basename(aln_fn).rpartition('_')[0]
    prefix = basename.partition('_')[0][:2]

    print basename
    
    aln = AlignIO.read(aln_fn, 'fasta')

    for i in range(aln.get_alignment_length()):
        print aln[:, i]
        counts = counts_from_col(aln[:, i])
        col_stats = column_stats(counts)
        
        print >>outfile,
        '\t'.join([ str(it) for it in [basename, i+1, col_stats[0], col_stats[1], entropy(counts)] ])
        
