"""
Download PDB files for all alignments. Choose the best structure for each
aligment and report it in the ensembl_pdb_map.txt.
"""
import os
import os.path as path
import glob
from argparse import ArgumentParser
from collections import defaultdict
import warnings
# Suppress warnings from PDBParser
warnings.simplefilter("ignore")

from Bio import SeqIO
from Bio import PDB

import utils
p = PDB.PDBParser()

argparser = ArgumentParser()

# argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--pdbdir', metavar='pdb_dir', type=str, required=True)
argparser.add_argument('--pdbmap', metavar='pdb_map', type=str, required=True)
argparser.add_argument('--log', metavar='log', type=str, required=True)
argparser.add_argument('--mismatch_log', metavar='mismatch_log', type=str, required=True)
argparser.add_argument('--missing_log', metavar='missing_log', type=str, required=True)

args = argparser.parse_args()


## Directory to store PDB files
pdb_dir = args.pdbdir
utils.check_dir(pdb_dir)

best_pdb_file = open(args.pdbmap)
logfile = open(args.log, 'w')
mismatch_file = open(args.mismatch_log, 'w')
missing_file = open(args.missing_log, 'w')


# We use the PDBe for faster transfers
# pdbl=PDB.PDBList(server="ftp://ftp.ebi.ac.uk/pub/databases/rcsb/")
pdbl=PDB.PDBList()
for l in best_pdb_file:
    ens, pdbname = l.rstrip().split('\t')[:2]
    final_file = None
    if path.exists(path.join(pdb_dir, 'pdb'+pdbname+'.ent')):
        print pdbname, "already downloaded."
        continue
    try:
        final_file = pdbl.retrieve_pdb_file(pdbname, pdir=pdb_dir)
    except IOError, e:
        if e.errno == "ftp error":
            print >>logfile, "Failed to download", pdbname
            continue
        else:
            print e.errno
            raise

    print pdbname, final_file
    basename = path.basename(final_file).partition('.')[0][-4:].lower()
    if basename != pdbname:
        print >>mismatch_file, pdbname, basename
