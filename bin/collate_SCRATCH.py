import sys
from os import path
import glob
import re
import argparse

import pandas
from Bio import SeqIO
import utils

argparser = argparse.ArgumentParser()

argparser.add_argument('--indir', metavar="input_dir", type=str, required=True)
argparser.add_argument('--outfile', metavar="output_dir", type=str, required=True)

def main():
    args = argparser.parse_args()

    ens_ids, seqsets, sites, buriedness, secstr = [], [], [], [], []

    files = glob.glob(path.join(args.indir, "*", "*.fa.acc"))
    print len(files)

    for f in files:
        dirname, basename = path.split(f)
        seqset = '_'.join(basename.split('.')[0].split('_')[1:])

        res = SeqIO.read(f, 'fasta')
        ss_f = SeqIO.read(open(path.join(dirname, path.splitext(basename)[0]+'.ss')), 'fasta')
        ens_id = res.id.partition(' ')[0]

        buried = [ 'exposed' if r == 'e' else 'buried' for r in res.seq ]
        ss = list(ss_f.seq)
        ens_ids.extend([ ens_id ]*len(buried))
        seqsets.extend([ seqset ]*len(buried))
        sites.extend(range(1, len(buried)+1))
        buriedness.extend(buried)
        secstr.extend(ss)
        
    out_df = pandas.DataFrame({ 'stable_id' : ens_ids,
                                'dataset' : seqsets,
                                'human_idx' : sites,
                                'buriedness_pred' : buriedness,
                                'ss_pred' : secstr })

    out_df.to_csv(args.outfile, index=False, sep="\t", quoting=False)

if __name__ == "__main__":
    main()
