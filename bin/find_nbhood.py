import re
import sys
import glob
from os import path
from argparse import ArgumentParser
import random
import utils
from collections import defaultdict
from pprint import pprint
from operator import itemgetter

from Bio import SeqIO
from Bio import AlignIO
from Bio import SeqUtils
from Bio import PDB

import pandas

re_resid = re.compile("(-?[0-9]+)([A-Z]*)")

def isnan(f):
    return not (f == f)

def parse_coord(coord):
    n, ic = re_resid.match(coord).groups()
    if ic == '':
        ic = ' '

    return ' ', int(n), ic


def compute_neighbours(start_res, chain, dist_thr):
    neighbours = set()

    for r1 in chain:
        if SeqUtils.seq1(r1.resname) == 'X':
            continue
        
        for at1 in r1:
            for r2 in start_res:
                if r1 == r2:
                    continue

                if SeqUtils.seq1(r2.resname) == 'X':
                    # Maybe kick up more of a fuss here if one of the residues of interest 
                    # is not a residue
                    continue
                    
                for at2 in r2:
                    if at1 - at2 < dist_thr:
                        neighbours.add(r1)

    return sorted(neighbours)
        
def parse_chunk(df, pdb_dir, dist):
    if len(set(df.pdb_chain)) > 1:
        print >>sys.stderr, "ERROR: multiple chains for a single ensembl ID"

    ens_id = df.stable_id.iloc[0]
    pdb_id = df.pdb_id.iloc[0]
    pdb_chain_id = df.pdb_chain.iloc[0]
    
    pdb_file = path.join(pdb_dir, 'pdb'+pdb_id+'.ent')
    pdb = p.get_structure(pdb_id, pdb_file)
    pdb_chain = pdb[0][pdb_chain_id]

    if len(pdb) > 1:
        print >>sys.stderr, "ERROR: multiple models in PDB file"

    print >>sys.stderr, pdb_file
    res_list = [ parse_coord(str(c)) for c in df.pdb_pos ]
    res_list = [ pdb_chain[c] for c in res_list if c in pdb_chain ]
    print >>sys.stderr, res_list
    nbhood = compute_neighbours(res_list, pdb_chain, dist)

    for res in [ ''.join([ str(c) for c in r.get_id() ]).strip() for r in nbhood ]:
        print '\t'.join([ ens_id, pdb_id, pdb_chain_id, res ])

    

p = PDB.PDBParser(QUIET=True)

argparser = ArgumentParser()
argparser.add_argument("--infile", metavar="input_file", type=str, required=True)
argparser.add_argument("--pdbdir", metavar="pdb_dir", type=str, required=True)
argparser.add_argument("--dist_thr", metavar="thr", type=float, default=4.0)

# TODO It would be better to read in a table with a header and then fish out the 
# columns of interest
def main():
    args = argparser.parse_args()
    
    infile = pandas.io.parsers.read_csv(args.infile, sep='\t')
    # infile.columns = [ 'ens_id', 'ens_pos', 'uni_id', 'uni_pos', 'pdb_id', 'pdb_chain', 'pdb_pos', 'ss', 'rsa', 'omega' ]

    infile.groupby(('stable_id', 'pdb_id')).apply(parse_chunk, args.pdbdir, args.dist_thr)

if __name__ == "__main__":
    main()
