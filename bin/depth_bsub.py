import os
from os import path
import subprocess
from argparse import ArgumentParser
import pandas

import utils

argparser = ArgumentParser()

argparser.add_argument('--infile', metavar='data_table', type=str, required=True)
argparser.add_argument('--pdbdir', metavar='pdb_root', type=str, required=True)
argparser.add_argument('--outdir', metavar='out_root', type=str, required=True)
argparser.add_argument('--logdir', metavar='log_file', type=str, required=True)
argparser.add_argument('--depth_cmd', metavar='depth_cmd', type=str, default="bin/run_depth.py")


if __name__ == "__main__":
    args = argparser.parse_args()

    infile = pandas.read_csv(open(args.infile), comment="\n", sep="\t")
    infile.drop_duplicates(subset=["pdb_id", "pdb_chain"], inplace=True)

    for i, r in infile.iterrows():
        p = subprocess.Popen(["bsub",
                              "-o"+path(args.logdir, outroot+'.log'),
                              "python " +
                              args.depth_cmd +
                              " --infile " + path.join(args.pdbdir, "pdb"+r['pdb_id']+".ent") +
                              " --stable_id " + r['stable_id'] +
                              " --pdb_id " + r['pdb_id'] +
                              " --pdb_chain " + r['pdb_chain'] +
                              " --outdir " + args.outdir +
                              " --logdir " + args.logdir])
        p.wait()
