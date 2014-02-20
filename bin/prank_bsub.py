from subprocess import Popen
from glob import glob
import os
from os import path
import argparse
import re

import utils

segfault_RE = re.compile("core dumped") # For rerunning

prank_cmd = "prank -d={} -t={} -o={} -prunetree -codon"

argparser = argparse.ArgumentParser()
argparser.add_argument('--clade', metavar='clade_name', type=str, required=True)
argparser.add_argument('--inroot', metavar='input_root', type=str, required=True)
argparser.add_argument('--treeroot', metavar='tree_root', type=str, required=True)
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)
argparser.add_argument('--rerun', action='store_true')

args = argparser.parse_args()

inroot = args.inroot
treeroot = args.treeroot
alndir = args.outroot
logroot = args.logdir

utils.check_dir(logroot)
utils.check_dir(alndir)

for infile in glob(path.join(inroot, clade, "*", "*.fa")):
    print infile
    basename = path.basename(infile).partition('.')[0]

    treedir = path.join(treeroot, clade, basename[:2])
    treefile = path.join(treedir, basename + '.nh')

    outdir = path.join(alndir, clade, basename[:2])
    utils.check_dir(outdir)
    outfile = path.join(outdir, basename + '_prank')

    logdir = path.join(logroot, clade, basename[:2])

    os.check_dir(logdir)

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
