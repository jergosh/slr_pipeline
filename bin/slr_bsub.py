from subprocess import Popen
from glob import glob
import os
from os import path
import argparse
import re

slr_cmd = "/nfs/research2/goldman/gregs/sw/slr-1.4.1/bin/Slr {}"

pr_root = "/nfs/research2/goldman/gregs/slr_pipeline"
inroot = path.join(pr_root, "data/ens/73/slr")
logroot = path.join(pr_root, "log/slr")

argparser = argparse.ArgumentParser()
args = argparser.parse_args()

clades = [ "Eutheria" ]
for clade in clades:
    for infile in glob(path.join(inroot, clade, "*", "*.ctl")):
    # for infile in glob(path.join(inroot, clade, "*", "*.ctl"))[:10]:
        basename = path.basename(infile)

        logdir = path.join(logroot, clade, basename[:2])

        if not path.exists(logdir):
            os.mkdir(logdir)

        logfile = path.join(logdir, basename.rpartition('.')[0] + '.log')
        # errfile = path.join(logdir, clade, basename + '.err')

        slr = slr_cmd.format(basename)

        p = Popen(['bsub', '-o'+logfile, '-g', '/slr',
                   '-cwd', path.dirname(infile), slr])
        
        p.wait()
