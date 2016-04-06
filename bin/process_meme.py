
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

species_RE = re.compile("([A-Z]+)")

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue)

argparser = argparse.ArgumentParser()

argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--memeroot', metavar='meme_root', type=str, required=True)
argparser.add_argument('--alnroot', metavar='aln_root', type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

args = argparser.parse_args()

clade = args.clade
alndir = args.alnroot
memeroot = args.memeroot
slr_all = args.outfile

all_ids = []
all_data = [] # pandas.DataFrame(columns=colnames)

for aln_fn in glob(path.join(alndir, clade, "*", "*_prank.best.fas")):
    dataset = path.basename(result_fn).partition('.')[0]
    prefix = basename.partition('_')[0][:2]

    result_fn = path.join(memeroot, clade, basename+'.txt')

    if not path.exists(result_fn):
        print result_fn, "doesn't exist!"
        continue

    result = pandas.read_table(open(result_fn))

    aln = AlignIO.read(aln_fn, 'fasta')
    for seqr in aln:
        species = species_RE.match(seqr.id).groups()[0]
        if species[:-1] != "ENS":
            continue

        all_ids.append(seqr.id)
        idx = [ i for (i, codon) in enumerate(grouper(seqr.seq, 3)) if ''.join(codon) != '---' ]
        result_subset = copy.deepcopy(result.ix[idx, :])
        result_subset.ix[:, 0] = idx
        result_subset.ix[:, 0] += 1

        meme_out = file(path.join(memeroot, clade, prefix, seqr.id + '_' + basename + '_meme.txt'), 'w')
        result_subset.to_csv(meme_out, quoting=False, index=False, sep='\t')

        result_subset['dataset'] = pandas.Series([dataset]*result_subset.shape[0], index=result_subset.index)
        result_subset['stable_id'] = pandas.Series([seqr.id]*result_subset.shape[0], index=result_subset.index)
        result_subset['human_idx'] = pandas.Series(range(1, result_subset.shape[0]+1), index=result_subset.index)

        all_data.append(result_subset)

    
all_data = pandas.concat(all_data)
all_data.to_csv(slr_all, quoting=False, index=False, sep='\t')
