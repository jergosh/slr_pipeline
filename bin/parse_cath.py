import sys
import pandas
import argparse
from os import path
from glob import glob
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from collections import defaultdict
from subprocess import Popen
import utils
from utils import parse_stk

argparser = argparse.ArgumentParser()

argparser.add_argument('--cath_file', metavar="file", type=str, required=True)
argparser.add_argument('--slr', metavar="file", type=str, required=True)
argparser.add_argument('--hmmdir', metavar="dir", type=str, required=True)
argparser.add_argument('--seqdir', metavar="dir", type=str, required=True)
argparser.add_argument('--outdir', metavar="dir", type=str, required=True)
argparser.add_argument('--outfile', metavar="file", type=str, required=True)

args = argparser.parse_args()

utils.check_dir(args.outdir)

print >>sys.stderr, "Loading SLR results...",
slr = pandas.read_table(open(args.slr), sep="\t")
print >>sys.stderr, "done."

id2ds = {}
def add_map(df):
    idx = df.index[0]
    id2ds[df.stable_id[idx]] = df.dataset[idx]
slr.groupby("stable_id").apply(add_map)

domain_dict = defaultdict(list)
for l in open(args.cath_file):
    f = l.rstrip().split('\t')
    
    domain_dict[f[2]].append(f)

all_ids = set(slr.stable_id)

outfile = open(args.outfile, 'w')
print >>outfile, '\t'.join(["stable_id", "dataset", "cath_id", "pdb_id", "pos", "omega"])

coords_dict = {}
for domain, domain_anns in domain_dict.items():
    out_seqs = []
    aln_fn = path.join(args.outdir, domain+'.fa')
    hmm_fn = path.join(args.hmmdir, domain_anns[0][3]+'.hmm')
    stockholm_fn = path.join(args.outdir, domain+'.stk')

    for domain_ann in domain_anns:
        ens, coords, cath_id, pdb_id = domain_ann 

        if ens not in all_ids:
            print "Skipping", ens
            continue
        else:
            print ens

        dataset = id2ds[ens]
        seq_fn = path.join(args.seqdir, dataset[:2], dataset+'.fa')
        seq = SeqIO.to_dict(SeqIO.parse(seq_fn, 'fasta'))[ens].seq

        domain_seq = Seq('')
        coords = [ int(i) for i in coords.split(':') ]
        coords_dict[(ens, cath_id)] = coords
        for begin, end in utils.grouper(coords, 2):
            domain_seq += seq[(begin-1):end] # CATH uses 1-based, inclusive coordinates

        if not len(domain_seq):
            continue

        out_seqs.append(SeqRecord(seq=domain_seq, id=ens, description=''))
    
    if not len(out_seqs):
        continue
    SeqIO.write(out_seqs, open(aln_fn, 'w'), 'fasta')

    hmmalign_cmd = ["hmmalign", "-o", stockholm_fn, hmm_fn, aln_fn ]

    p = Popen(hmmalign_cmd)
    p.wait()

    # try:
    hmm_aln = parse_stk(open(stockholm_fn)) # AlignIO.read(stockholm_fn, 'stockholm')
    print hmm_aln
    # except ValueError, e:
    #     continue
    
    for seq_id, seq in hmm_aln.items():
        seq_id = seq_id.partition('_')[0]
        slr_subset = slr[slr.stable_id == seq_id]
        slr_subset.reset_index(inplace=True)

        omegas_all = slr_subset.Omega
        omegas = []
        dataset = id2ds[seq_id]
        i_seq = 0

        print seq_id, coords_dict[(seq_id, cath_id)]
        print coords_dict[(seq_id, cath_id)]
        for begin, end in utils.grouper(coords_dict[(seq_id, cath_id)], 2):
            omegas.extend(list(omegas_all[(begin-1):end])) # CATH uses 1-based, inclusive coordinates
        print len(omegas)
        for i, c in enumerate(seq):
            i = str(i+1)
            if c.isupper():
                omega = str(omegas[i_seq])
                i_seq += 1
            elif c == "-":
                omega = 'NA'
            elif c.islower():
                i_seq += 1
                continue
            elif c == ".":
                continue
            else:
                print "ERROR"
                sys.exit(-1)

            print >>outfile, '\t'.join([ seq_id, dataset, cath_id, pdb_id, i, omega ])

        print i_seq, len(omegas)
        # assert i_seq == len(omegas)
