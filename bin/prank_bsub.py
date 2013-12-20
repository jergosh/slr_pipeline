from subprocess import Popen
from glob import glob
import os
from os import path

prank_cmd = "prank -d={} -t={} -o={} -prunetree"

pr_root = "/nfs/research2/goldman/gregs/slr_pipeline"
indir = path.join(pr_root, "data/ens/73/tmp")
treedir = path.join(pr_root, "data/ens/73/seqsets")
alndir = path.join(pr_root, "data/ens/73/aln")
logdir = path.join(pr_root, "log/prank")

clades = [ "Eutheria" ]
for clade in clades:
    for infile in glob(path.join(indir, clade, "*.fa")):
    # for infile in glob(path.join(indir, clade, "*.fa"))[:2]:
        basename = path.basename(infile).partition('.')[0]

        treefile = path.join(treedir, clade, basename + '.nh')

        outdir = path.join(alndir, clade, basename[:2])
        if not os.exists(outdir):
            os.mkdir(outdir)
        outfile = path.join(outdir, basename + '_prank')

        if not os.exists(logdir):
            os.mkdir(logdir)
        logdir = path.join(logdir, clade, basename[:2])

        logfile = path.join(logdir, basename + '.log')
        errfile = path.join(logdir, clade, basename + '.err')

        prank = prank_cmd.format(infile, treefile, outfile)
        p = Popen(['bsub', '-R "rusage[tmp=512]"' '-o'+logfile, '-cwd /tmp',
                   '-g /prank', prank])
        p.wait()
