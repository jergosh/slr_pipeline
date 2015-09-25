from subprocess import Popen
from glob import glob
import os
from os import path
import re

segfault_RE = re.compile("core dumped")

prank_cmd = "prank -d={} -t={} -o={} -prunetree"

pr_root = "/nfs/research2/goldman/gregs/slr_pipeline"
indir = path.join(pr_root, "data/ens/73/tmp")
treedir = path.join(pr_root, "data/ens/73/seqsets")
alndir = path.join(pr_root, "data/ens/73/aln")
logdir = path.join(pr_root, "log/prank")

clades = [ "Eutheria" ]
for clade in clades:
    for logfile in glob(path.join(logdir, clade, "*.log")):
        basename = path.basename(logfile).partition('.')[0]
        if segfault_RE.search(open(logfile).read()) != None:
            treefile = path.join(treedir, clade, basename + '.nh')
            infile = path.join(indir, clade, basename + '.fa')
            outfile = path.join(alndir, clade, basename + '_prank')
            errfile = path.join(logdir, clade, basename + '.err')
            os.remove(logfile)

            prank = prank_cmd.format(infile, treefile, outfile)
            p = Popen(['bsub', '-o'+logfile, '-cwd /tmp', prank])
            p.wait()
