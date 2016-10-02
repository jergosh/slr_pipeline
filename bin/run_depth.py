import os
from os import path
import subprocess
from argparse import ArgumentParser

import utils

argparser = ArgumentParser()

argparser.add_argument('--infile', metavar='pdb_file', type=str, required=True)
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
argparser.add_argument('--logfile', metavar='log_file', type=str, required=True)
argparser.add_argument('--depth_exec', metavar='bin', type=str,
                       default="/hps/nobackup/goldman/gregs/depth-2.0/doc")


if __name__ == "__main__":
    args = argparser.parse_args()

    depth_cmd = "{} -i {} -o {}".format(args.depth_exec, args.infile, args.outroot)
    
    p = subprocess.Popen[ "bsub", "-o"+args.logfile, depth_cmd ]
    p.wait()

    for l in open(args.outroot + "-residue.depth"):
        
