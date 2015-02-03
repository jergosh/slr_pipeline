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
yeast_RE = re.compile("Y[A-P][LR][0-9]{3}[WC]")

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
slr_all = args.outfile

all_ids = []
all_data = pandas.DataFrame(columns=colnames)

for aln_fn in glob(path.join(alndir, clade, "*", "*_prank.best.fas")):
    basename = path.basename(aln_fn).rpartition('_')[0]
    prefix = basename.partition('_')[0][:2]

    slr_fn = path.join(slrroot, clade, prefix, basename+'_matched.res')

    if not path.exists(slr_fn):
        print slr_fn, "doesn't exist!"
        continue

    # TODO Make sure colspecs work in all cases
    slr = pandas.read_fwf(open(slr_fn), colspecs=colspecs, header=False, comment="\n")

    aln = AlignIO.read(aln_fn, 'fasta')
    for seqr in aln:
        if args.clade == "yeast":
            if yeast_RE.match(seqr.id) is None:
                continue
        else:
            species = species_RE.match(seqr.id).groups()[0]
            if species[:-1] != "ENS":
                continue

        all_ids.append(seqr.id)
        idx = [ i for (i, codon) in enumerate(grouper(seqr.seq, 3)) if ''.join(codon) != '---' ]
        slr_subset = copy.deepcopy(slr.ix[idx, :])
        slr_subset.ix[:, 0] = idx
        slr_subset.ix[:, 0] += 1

        slr_out = file(path.join(slrroot, clade, prefix, seqr.id + '_' + basename + '_matched.res'), 'w')
        slr_subset.to_csv(slr_out, quoting=False, index=False, sep='\t')

        # slr_subset.insert(0, 'dataset', pandas.Series([basename]*slr_subset.shape[0]))
        # slr_subset.insert(0, 'stable_id', pandas.Series([seqr.id]*slr_subset.shape[0]))
        slr_subset['dataset'] = pandas.Series([basename]*slr.shape[0])
        slr_subset['stable_id'] = pandas.Series([seqr.id]*slr.shape[0])
        slr_subset['human_idx'] = pandas.Series(range(1, slr.shape[0]+1))

        all_data = all_data.append(slr_subset)

all_data.rename(columns={"# Site": "Site"}, inplace=True)
all_data.to_csv(slr_all, quoting=False, index=False, sep='\t')
