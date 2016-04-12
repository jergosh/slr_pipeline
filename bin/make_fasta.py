import sys
import re 
import os
from os import path
from glob import glob
import pickle
from argparse import ArgumentParser

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import utils

import ete3

transcript_RE = re.compile("transcript:(\w+)")

argparser = ArgumentParser()

argparser.add_argument('--cds', metavar='cds_dir', type=str, required=True)
argparser.add_argument('--pep', metavar='pep_dir', type=str, required=True)
argparser.add_argument('--infile', metavar='inroot', type=str, required=True)
argparser.add_argument('--outfile', metavar='outroot', type=str, required=True)

args = argparser.parse_args()

print >>sys.stderr, "Building the ID dictionary..."
t2p = {} # Transcript to protein ID
ens_pep_dir = args.pep
for f in glob(path.join(ens_pep_dir, '*.all.fa')):
    print "Processing", f
    for seqr in SeqIO.parse(f, 'fasta'):
        tr = transcript_RE.search(seqr.description)
        if not tr:
            print seqr.description
            raise ValueError("Transcript ID not found in FASTA description field.")

        tr_name = tr.groups()[0]
        t2p[tr_name] = seqr.id

## Build a mapping of Ens IDs to coding sequences 
ens_map = {}
ens_cdna_dir = args.cds
for f in glob(path.join(ens_cdna_dir, '*.fa')):
    print "Processing", f
    for seqr in SeqIO.parse(f, 'fasta'):
        if seqr.id in ens_map:
            print "Duplicate id", seqr.id
            sys.exit(-1)

        ens_map[t2p[seqr.id]] = seqr.seq

tree = ete3.Tree(args.infile)

seqs = []
for l in tree.iter_leaves():
    seq = ens_map.get(seqid)
    if seq is None:
        print >>sys.stderr, seqid, "missing!"
    else:
        seqs.append(SeqRecord(seq, id=seqid, description=""))

SeqIO.write(seqs, open(args.outfile, 'w'), 'fasta')
