import os
from os import path
import subprocess
from argparse import ArgumentParser

import utils

argparser = ArgumentParser()

argparser.add_argument('--infile', metavar='pdb_file', type=str, required=True)
argparser.add_argument('--stable_id', metavar='stable_id', type=str, required=True)
argparser.add_argument('--pdb_id', metavar='pdb_id', type=str, required=True)
argparser.add_argument('--pdb_chain', metavar='pdb_chain', type=str, required=True)
argparser.add_argument('--outdir', metavar='out_dir', type=str, required=True)
argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)
argparser.add_argument('--depth_exec', metavar='bin', type=str,
                       default="/hps/nobackup/goldman/gregs/depth-2.0/bin/DEPTH")


if __name__ == "__main__":
    args = argparser.parse_args()
    outroot = "_".join([args.stable_id, args.pdb_id, args.pdb_chain])

    depth_cmd = "{} -i {} -o {}".format(args.depth_exec, args.infile,
                                        path.join(args.outdir, args.outroot))
    
    p = subprocess.Popen([ "bsub", "-o"+path(args.logdir, outroot+'.log'), depth_cmd ])
    p.wait()

    outfile = open(path.join(args.outdir,
                             outroot+'.tab'), 'w')
    print >>outfile, "\t".join(["stable_id", "pdb_id", "pdb_chain", "pdb_pos", "depth"])

    infile = open(path.join(args.outdir, outroot) + "-residue.depth")
    infile.readline()
    for l in infile:
        f = l.rstrip().split('\t')
        pdb_chain, pdb_pos = f[0].split(':')

        if pdb_chain == args.pdb_chain:
            print >>outfile, "\t".join([args.stable_id, args.pdb_id, pdb_chain, pdb_pos, f[2]])
