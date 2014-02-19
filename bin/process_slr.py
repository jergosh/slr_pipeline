from glob import glob
import os
from os import path
import itertools
import re
from Bio import AlignIO 
import pandas
import sys
import copy

species_RE = re.compile("([A-Z]+)")

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue)

pr_root = "/nfs/research2/goldman/gregs/slr_pipeline"
inroot = path.join(pr_root, "data/ens/73/seqsets_cds")
treeroot = path.join(pr_root, "data/ens/73/seqsets")
alndir = path.join(pr_root, "data/ens/73/aln")
slrroot = path.join(pr_root, "data/ens/73/slr")
slr_all = path.join(pr_root, "data/ens/73/slr_all.tab")

logroot = path.join(pr_root, "log/prank")

colspecs = [(0, 8), (8, 17), (17, 25), (25, 34), (34, 43), (43, 51), (51, 60), (60, 71), (71, 82), (82, 92), (92, 99), (99, 111)]
colnames = [ "dataset", "stable_id", "human_idx", "# Site", "Neutral", "Optimal", "Omega", "lower", "upper", "LRT_Stat", "Pval",
             "Adj.Pval", "Q-value", "Result", "Note" ]
all_ids = []
all_data = pandas.DataFrame(columns=colnames)

clade = "Eutheria"
for aln_fn in glob(path.join(alndir, clade, "*", "*_prank.best.fas")):
    basename = path.basename(aln_fn).rpartition('_')[0]

    slr_fn = path.join(slrroot, clade, basename[:2], basename+'_matched.res')

    if not path.exists(slr_fn):
        print slr_fn, "doesn't exist!"
        continue

    slr = pandas.read_fwf(open(slr_fn), colspecs=colspecs, header=False, comment="\n")

    aln = AlignIO.read(aln_fn, 'fasta')
    for seqr in aln:
        species = species_RE.match(seqr.id).groups()[0]
        if species[:-1] != "ENS":
            continue

        all_ids.append(seqr.id)
        idx = [ i for (i, codon) in enumerate(grouper(seqr.seq, 3)) if ''.join(codon) != '---' ]
        slr_subset = copy.deepcopy(slr.ix[idx, :])
        slr_subset.ix[:, 0] = idx
        slr_subset.ix[:, 0] += 1

        slr_out = file(path.join(slrroot, clade, basename[:2], seqr.id + '_' + basename + '_matched.res'), 'w')
        slr_subset.to_csv(slr_out, quoting=False, index=False, sep='\t')

        # slr_subset.insert(0, 'dataset', pandas.Series([basename]*slr_subset.shape[0]))
        # slr_subset.insert(0, 'stable_id', pandas.Series([seqr.id]*slr_subset.shape[0]))
        slr_subset['dataset'] = pandas.Series([basename]*slr.shape[0])
        slr_subset['stable_id'] = pandas.Series([seqr.id]*slr.shape[0])
        slr_subset['human_idx'] = pandas.Series(range(1, slr.shape[0]+1))

        all_data = all_data.append(slr_subset)

all_data.rename(columns={"# Site": "Site"}, inplace=True)
all_data.to_csv(slr_all, quoting=False, index=False, sep='\t')
