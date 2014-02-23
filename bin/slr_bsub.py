from subprocess import Popen
from glob import glob
import os
from os import path
import argparse
import re
import utils

slr_cmd = "/nfs/research2/goldman/gregs/sw/slr-1.4.1/bin/Slr {}"

argparser = argparse.ArgumentParser()

argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--slrroot', metavar='slr_root', type=str, required=True)
argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)

args = argparser.parse_args()

inroot = args.slrroot
logroot = args.logdir

utils.check_dir(logroot)

# for infile in glob(path.join(inroot, args.clade, "*", "*.ctl")):
for infile in glob(path.join(inroot, args.clade, "*", "*.ctl"))[:10]:
    basename = path.basename(infile)

    logdir = path.join(logroot, args.clade, basename[:2])

    utils.check_dir(logdir)

    logfile = path.join(logdir, basename.rpartition('.')[0] + '.log')
    # errfile = path.join(logdir, args.clade, basename + '.err')

    slr = slr_cmd.format(basename)

    p = Popen(['bsub', '-o'+logfile, '-g', '/slr',
               '-cwd', path.dirname(infile), slr])
        
    p.wait()
