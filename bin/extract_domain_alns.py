import glob
import argparse
from os import path
import pandas

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_protein, Gapped
import utils


argparser = argparse.ArgumentParser()

argparser.add_argument('--indir', metavar="dir", type=str, required=True)
argparser.add_argument('--domaindir', metavar="dir", type=str, required=True)
argparser.add_argument('--dataset_map', metavar="file", type=str, required=True)
argparser.add_argument('--cath', metavar="file", type=str, required=True)
argparser.add_argument('--outdir', metavar="output_dir", type=str, required=True)

def extract_aln(df, aln_dir, domain_dir, dataset_map, outdir):
    cath_id = df.cath_id.iloc[0]
    stable_id = df.stable_id.iloc[0]

    print cath_id, stable_id
    if stable_id not in dataset_map:
        return

    dataset = dataset_map[stable_id]
    subdir = dataset.split('_')[0][:2]

    aln_file = path.join(aln_dir, subdir, dataset+'_prank.best.fas')
    aln_dict = SeqIO.to_dict(SeqIO.parse(aln_file, 'fasta'))
    ens_seq = aln_dict[stable_id]
        
    human_seqs = []
    for seqr in aln_dict.values():
        out_seq = []
        for i, c in enumerate(ens_seq):
            if c != '-':
                out_seq.append(seqr.seq[i])

        out_seq = ''.join(out_seq)
        human_seqs.append(SeqRecord(Seq(out_seq), id=seqr.id, description=seqr.description))

    stk_aln_file = path.join(domain_dir, cath_id+'.stk')
    stk_aln = utils.parse_stk(open(stk_aln_file))
    domain_seq = stk_aln[stable_id]
    i_seq = 0
    out_seqs = []
    for seqr in human_seqs:
        out_seqs.append(SeqRecord(Seq(""), id=seqr.id, description=seqr.description))

    out_seqs = MultipleSeqAlignment(out_seqs) #, Gapped(generic_protein, '-'))
    human_seqs = MultipleSeqAlignment(human_seqs)

    for i, c in enumerate(domain_seq):
        if c.isupper():
            out_seqs = out_seqs + human_seqs[:, i_seq:i_seq+1]
            i_seq += 1
        elif c == "-":
            for seqr in out_seqs.get_all_seqs():
                seqr.seq += '-'
        elif c.islower():
            i_seq += 1
            continue
        elif c == ".":
            continue
        else:
            print "ERROR"
            sys.exit(-1)
        
    outdir = path.join(outdir, subdir)
    outfile = path.join(outdir, stable_id+'_'+dataset+'.fa')
    utils.check_dir(outdir)
    SeqIO.write(out_seqs, open(outfile, 'w'), 'fasta')

    

def main():
    args = argparser.parse_args()

    utils.check_dir(args.outdir)

    dataset_map = {}
    for l in open(args.dataset_map):
        f = l.rstrip().split('\t')
        dataset_map[f[0]] = f[1]

    cath = pandas.read_table(args.cath, sep='\t')
    cath.groupby(['cath_id', 'stable_id']).apply(extract_aln, args.indir, args.domaindir, dataset_map, args.outdir)


if __name__ == "__main__":
    main()
    
