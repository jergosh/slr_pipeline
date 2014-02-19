from glob import glob
import os
from os import path
import re
import sys
import numpy as np
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pr_root = "/nfs/research2/goldman/gregs/slr_pipeline"
inroot = path.join(pr_root, "data/ens/73/seqsets_cds")
treeroot = path.join(pr_root, "data/ens/73/seqsets")

n_human = []
aln_sizes = []
paralog_fraction = []
clade = "Eutheria"
for treefile in glob(path.join(treeroot, clade, '*', '*.nh')):
    basename = path.basename(treefile).rpartition('.')[0]
    dirname = path.dirname(treefile)
    subsetfile = path.join(dirname, basename + '.tab')
    
    species = defaultdict(int)
    n = 0
    for l in open(subsetfile):
        n += 1
        f = l.rstrip().split('\t')
        species[f[1]] += 1

    paralog_fraction.append(float(n-len(species))/n)
    aln_sizes.append(n)
    n_human.append(species["homo sapiens"])

pages = PdfPages("subsets_stats.pdf")
fig = plt.figure()
n, bins, patches = plt.hist(n_human, 50, normed=0, histtype='stepfilled')
plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
plt.title("Number of human genes in alignment")
plt.xlabel("# of human genes")
plt.ylabel("Count")
pages.savefig(fig, bbox_inches=0)

fig = plt.figure()
n, bins, patches = plt.hist(aln_sizes, 50, normed=0, histtype='stepfilled')
plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
plt.title("Number of sequences in alignment")
plt.xlabel("# of sequences")
plt.ylabel("Count")
pages.savefig(fig, bbox_inches=0)

fig = plt.figure()
n, bins, patches = plt.hist(paralog_fraction, 50, normed=0, histtype='stepfilled')
plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
plt.title("Fraction of paralogs in aligment")
plt.xlabel("Fraction of paralogs")
plt.ylabel("Count")
pages.savefig(fig, bbox_inches=0)

pages.close()
