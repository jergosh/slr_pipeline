import re
import sys
import glob
import csv
from os import path
from argparse import ArgumentParser
from subprocess import Popen
import operator

from Bio import PDB
import numpy as np
import pandas

def rid2str(r):
    return r.id[0] + str(r.id[1]) + r.id[2]


argparser = ArgumentParser()
argparser.add_argument("--pdbmap", metavar="pdb_map", type=str, required=True)
argparser.add_argument("--pdbdir", metavar="pdb_dir", type=str, required=True)
argparser.add_argument("--outfile", metavar="out_dir", type=str, required=True)
argparser.add_argument("--thr", metavar="thr", type=float, default=0.05)

args = argparser.parse_args()

def find_neighbourhood(df, thr):
    pdb_id = df.pdb_id.iloc[0]
    stable_id = df.stable_id.iloc[0]
    chain_id = df.pdb_chain.iloc[0]

    try:
        pdb = p.get_structure(pdb_id, pdbfile)
        pdb_chain = pdb[0][chain_id]
    except IOError, e:
        print >>sys.stderr, "PDB file", pdb_id, "missing!"
        return

    out_df = pandas.DataFrame(columns=('source', 'pdb_pos'))
    out_residues = set()
    for i, row in df.iterrows():
        res_id = parse_coord(row.pdb_pos)
        try:
            r = pdb_chain[res_id]
        except KeyError, e:
            r = find_sequential(pdb_chain, res_id)
            if r is None:
                raise e

        for r2 in pdb_chain:
            if r2 in out_residues:
                print "Residue already in output_set"
                continue
            
            if r - r2 < thr:
                out_df.loc[df.shape[0]] = [ rid2str(r), rid2str(r2) ]
                out_residues.add(r2)

    return out_df

if __name__ == "__main__":
    pdb_map = pandas.read_table(args.pdbmap, dtype={ "stable_id": str, "pdb_id": str, "pdb_pos": str, "omega": np.float64 })

    out_map = pdb_map.groupby(["stable_id", "pdb_id", "pdb_chain"]).apply(find_neighbourhood, args.thr)
    out_map.to_csv(args.outfile, sep="\t", quoting=csv.QUOTE_NONE)

