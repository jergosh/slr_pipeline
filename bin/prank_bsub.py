from subprocess import Popen
from glob import glob
import os
from os import path
import argparse
import re

segfault_RE = re.compile("core dumped") # For rerunning

prank_cmd = "prank -d={} -t={} -o={} -prunetree -codon"

pr_root = "/nfs/research2/goldman/gregs/slr_pipeline"
inroot = path.join(pr_root, "data/ens/73/seqsets_cds")
treeroot = path.join(pr_root, "data/ens/73/seqsets")
alndir = path.join(pr_root, "data/ens/73/aln")
logroot = path.join(pr_root, "log/prank")

argparser = argparse.ArgumentParser()
argparser.add_argument('--rerun', action='store_true')

args = argparser.parse_args()

clades = [ "Eutheria" ]
for clade in clades:
    for infile in glob(path.join(inroot, clade, "*", "*.fa")):
        print infile
        basename = path.basename(infile).partition('.')[0]

        treedir = path.join(treeroot, clade, basename[:2])
        treefile = path.join(treedir, basename + '.nh')

        outdir = path.join(alndir, clade, basename[:2])
        if not path.exists(outdir):
            os.mkdir(outdir)
        outfile = path.join(outdir, basename + '_prank')

        logdir = path.join(logroot, clade, basename[:2])

        if not path.exists(logdir):
            os.mkdir(logdir)

        logfile = path.join(logdir, basename + '.log')
        errfile = path.join(logdir, clade, basename + '.err')

        prank = prank_cmd.format(infile, treefile, outfile)

        if args.rerun:
            if segfault_RE.search(open(logfile).read()) is None:
                continue
            else:
                os.remove(logfile)

        p = Popen(['bsub', '-R', 'rusage[tmp=512]', '-o'+logfile, '-g', '/prank',
                   '-cwd /tmp', prank])
        
        p.wait()
