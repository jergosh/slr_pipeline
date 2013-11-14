import sys
import re 
from os import path
from glob import glob
import pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

transcript_RE = re.compile("transcript:(\w+)")

pr_root = "."

all_species = []
species_file = "data/specieslist.txt"
for l in open(species_file):
    f = l.strip().split('\t')
    all_species.append(f[1].lower())
    print all_species[-1]


print >>sys.stderr, "Building the ID dictionary..."
## Build a mapping of Ens IDs to protein sequences 
t2p = {} # Transcript to protein ID
ens_pep_dir = path.join(pr_root, "data/ens/73/pep")
for f in glob(path.join(ens_pep_dir, '*.all.fa')):
    for seqr in SeqIO.parse(f, 'fasta'):
        tr = transcript_RE.search(seqr.description)
        if not tr:
            print seqr.description
            raise ValueError("Transcript ID not found in FASTA description field.")

        tr_name = tr.groups()[0]
        t2p[tr_name] = seqr.id

ens_map = {}
ens_cdna_dir = path.join(pr_root, "data/ens/73/cds/all")
for f in glob(path.join(ens_cdna_dir, '*_cds.fa')):
    for seqr in SeqIO.parse(f, 'fasta'):
        if seqr.id in t2p:
            pid = t2p[seqr.id]

            if pid in ens_map:
                print "Duplicate id", pid
                sys.exit(-1)

            ens_map[pid] = seqr.seq

clades_pickle = "data/clades.pk"
clades = pickle.load(open(clades_pickle))

indir = "data/ens/73/seqsets"
outdir = "data/ens/73/tmp"
clades = [ "Eutheria" ]



for clade in clades:
    for seqset in glob(path.join(indir, clade, "*.tab")):
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
            SeqIO.write(seqs, open(path.join(outdir, clade, setid + '.fa'), 'w'), 'fasta')