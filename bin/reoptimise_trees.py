import os 
from os import path as path
import shutil
import glob
import re
import utils
import argparse
from subprocess import Popen

raxml_cmd = "src/standard-RAxML-8.2.4/raxmlHPC-PTHREADS-SSE3 -T 1 -f e -m PROTGAMMAWAG -n {} \
-s {} -t {} -w {}"

argparser = argparse.ArgumentParser()
argparser.add_argument('--alndir', metavar='tree_directory', type=str, required=True)
argparser.add_argument('--treedir', metavar='tree_directory', type=str, required=True)
argparser.add_argument('--outroot', metavar='output_directory', type=str, required=True)
argparser.add_argument('--logroot', metavar='log_directory', type=str, required=True)
args = argparser.parse_args()

utils.check_dir(args.outroot)
utils.check_dir(args.logroot)

for aln_file in glob.glob(path.join(args.alndir, "*.fa")):
    print aln_file

    basename = path.basename(aln_file)
    input_name = basename.rpartition('_')[0]
    prefix = input_name.partition('_')[0][:2]
    tree_file = path.join(args.treedir, prefix, input_name+'.nh')

    logdir = path.join(args.logroot, prefix)
    utils.check_dir(logdir)
    logfile = path.join(logdir, input_name+'.log')

    outdir = path.join(args.outroot, prefix)
    utils.check_dir(outdir)
    raxml = raxml_cmd.format(input_name, path.abspath(aln_file),
                             path.abspath(tree_file), path.abspath(outdir))

    p = Popen([ 'bsub', '-o'+logfile, raxml ])
    p.wait()
