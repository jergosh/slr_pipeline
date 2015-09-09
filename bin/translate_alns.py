import glob
import argparse
from os import path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Data import CodonTable

import utils


argparser = argparse.ArgumentParser()

argparser.add_argument('--indir', metavar="dir", type=str, required=True)
argparser.add_argument('--outdir', metavar="dir", type=str, required=True)

def main():
    args = argparser.parse_args()

    codontable = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
    
    utils.check_dir(args.outdir)
    for infile in glob.glob(path.join(args.indir, '*', '*.best.fas')):
        dir, basename = path.split(infile)
        _, subdir = path.split(dir)
        print basename

        outseqs = []
        for seq in SeqIO.parse(open(infile), 'fasta'):
            outstr = ''.join([ codontable[''.join(c)] if ''.join(c) not in [ '---', 'NNN' ] else '-' 
                               for c in utils.grouper(seq, 3) ])
            outseq = SeqRecord(Seq(outstr), id=seq.id, description=seq.description)
            outseqs.append(outseq)

        outdir = path.join(args.outdir, subdir)
        utils.check_dir(outdir)
        SeqIO.write(outseqs, open(path.join(outdir, basename), 'w'), 'fasta')

if __name__ == "__main__":
    main()
