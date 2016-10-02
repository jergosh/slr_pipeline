import os
from os import path
import subprocess
from argparse import ArgumentParser
import pandas

import utils

argparser = ArgumentParser()

argparser.add_argument('--infile', metavar='data_table', type=str, required=True)
argparser.add_argument('--outdir', metavar='out_root', type=str, required=True)
argparser.add_argument('--logdir', metavar='log_file', type=str, required=True)


if __name__ == "__main__":
    args = argparser.parse_args()

    infile = pandas.read_csv(open(args.infile), comment="\n", sep="\t")
    infile.drop_duplicates(cols=["pdb_id", "pdb_chain"], inplace=True)
    print infile
