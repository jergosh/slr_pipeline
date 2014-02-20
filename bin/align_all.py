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

argparser = ArgumentParser()

argparser.add_argument('--cds', metavar='cds_dir', type=str, required=True)
argparser.add_argument('--pep', metavar='pep_dir', type=str, required=True)
argparser.add_argument('--inroot', metavar='inroot', type=str, required=True)
argparser.add_argument('--outroot', metavar='outroot', type=str, required=True)
argparser.add_argument('--specieslist', metavar='species_file', type=str, required=True)

args = argparser.parse_args()

transcript_RE = re.compile("transcript:(\w+)")

all_species = []
species_file = args.specieslist
for l in open(species_file):
    f = l.strip().split('\t')
    all_species.append(f[1].lower())
    print all_species[-1]


print >>sys.stderr, "Building the ID dictionary..."
t2p = {} # Transcript to protein ID
ens_pep_dir = path.join(pr_root, "data/ens/73/pep")
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
ens_cdna_dir = path.join(pr_root, "data/ens/73/cds/canonical")
for f in glob(path.join(ens_cdna_dir, '*_cds.fa')):
    print "Processing", f
    for seqr in SeqIO.parse(f, 'fasta'):
        if seqr.id in ens_map:
            print "Duplicate id", seqr.id
            sys.exit(-1)

        ens_map[t2p[seqr.id]] = seqr.seq

clades_pickle = "data/clades.pk"
clades = pickle.load(open(clades_pickle))

inroot = args.inroot
outroot = args.outroot


for seqset in glob(path.join(inroot, args.clade, "*", "*.tab")):
    setid = path.basename(seqset).rpartition('.')[0]
    seqs = []
    for l in open(seqset):
        seqid, species = l.rstrip().split('\t')
        if species not in all_species:
            continue

        seq = ens_map.get(seqid)
        if seq is None:
            print seqid, "missing!"
            break
        seqs.append(SeqRecord(seq, id=seqid, description=""))
    else:
        outdir =  path.join(outroot, args.clade, setid[:2])
        if not path.exists(outdir):
            os.mkdir(outdir)

        SeqIO.write(seqs, open(path.join(outdir, setid + '.fa'), 'w'), 'fasta')
