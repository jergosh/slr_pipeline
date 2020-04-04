from glob import glob
import os
from os import path
import itertools
import re
from Bio import AlignIO 
import pandas
import sys
import copy
import argparse

from slr import *

species_RE = re.compile("([A-Z]+)")

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue)

argparser = argparse.ArgumentParser()

argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--slrroot', metavar='slr_root', type=str, required=True)
argparser.add_argument('--alnroot', metavar='aln_root', type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

args = argparser.parse_args()

clade = args.clade
alndir = args.alnroot
slrroot = args.slrroot

all_ids = []
all_data = [] # pandas.DataFrame(columns=colnames)

outfile = open(args.outfile, 'w')
for aln_fn in glob(path.join(alndir, clade, "*", "*_prank.best.fas")):
    basename = path.basename(aln_fn).rpartition('_')[0]
    prefix = basename.partition('_')[0][:2]

    slr_fn = path.join(slrroot, clade, prefix, basename+'_matched.res')
    if not path.exists(slr_fn):
        print slr_fn, "doesn't exist!"
        continue

    slr = open(slr_fn)
    slr.readline()
    print >>outfile, slr.readline()
