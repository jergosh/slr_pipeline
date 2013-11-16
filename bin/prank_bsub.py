from subprocess import Popen
from glob import glob
from os import path

prank_cmd = "prank -d={} -t={} -o={} -prunetree"

pr_root = "/nfs/research2/goldman/gregs/slr_pipeline"
indir = path.join(pr_root, "data/ens/73/tmp")
treedir = path.join(pr_root, "data/ens/73/seqsets")
alndir = path.join("data/ens/73/aln")
logdir = path.join(pr_root, "log/prank")

clades = [ "Eutheria" ]
for clade in clades:
    for infile in glob(path.join(indir, clade, "*.fa")):
    # for infile in glob(path.join(indir, clade, "*.fa"))[:2]:
        basename = path.basename(infile).partition('.')[0]

        treefile = path.join(treedir, clade, basename + '.nh')
        outfile = path.join(alndir, clade, basename + '_prank')
        logfile = path.join(logdir, clade, basename + '.log')
        errfile = path.join(logdir, clade, basename + '.err')

        prank = prank_cmd.format(infile, treefile, outfile)
        p = Popen(['bsub', '-o'+logfile, '-cwd /tmp', prank])
        p.wait()
