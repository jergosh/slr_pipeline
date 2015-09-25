import os
from os import path
import sys
import argparse
import pandas
import utils

argparser = argparse.ArgumentParser()

argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--infile', metavar='input_file', type=str, required=True)
argparser.add_argument('--outdir', metavar='output_root', type=str, required=True)

args = argparser.parse_args()

utils.check_dir(args.outdir)
outdir = path.join(args.outdir, args.clade)
utils.check_dir(outdir)

def write_df(df):
    outfile = path.join(outdir, df.reset_index().loc[0, "ens"]+'.tab')
    print outfile
    df.to_csv(outfile, quoting=False, index=False, sep='\t')
    
pdb_master_table = pandas.io.parsers.read_csv(args.infile, sep='\t', header=False)
pdb_master_table.columns = ['ens', 'pdb', 'pos', 'secondary', 'rsa', 'omega']
pdb_master_table.groupby('ens').apply(write_df)

