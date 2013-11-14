from subprocess import Popen
from glob import glob
from os import path

prank_cmd = "prank -d={} -t={} -o={} -prunetree"

indir = "data/ens/73/tmp/"
treedir = "data/ens/73/seqsets"
alndir = "data/ens/73/aln"

# for infile in glob(path.join(indir, "*.fa")):
for infile in glob(path.join(indir, "*.fa"))[:2]:
    basename = path.basename(infile).partition('.')[0]

    treefile = path.join(treedir, basename + '.nh')
    outfile = path.join(alndir, basename + '_prank')

    prank = prank_cmd.format(infile, treefile, outfile)

    p = Popen('bsub', '-o /deb/null', prank)
    p.wait()
