import subprocess
import os.path as path
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Pass the genome file as reference and a list of IDs?

argparser = argparse.ArgumentParser()

argparser.add_argument('--infile', metavar="file", type=str, required=True)
argparser.add_argument('--reffile', metavar="file", type=str, required=True)
argparser.add_argument('--outdir', metavar="dir", type=str, required=True)
argparser.add_argument('--iupre', metavar="iupred", type=str, default="iupred")


if __name__ == "__main__":
    args = argparser.parse_args()

    ref_handle = open(args.reffile)
    ref_dict = SeqIO.to_dict(SeqIO.parse(ref_handle, "fasta"))

    infile = open(args.infile)

    for l in infile:
        f = l.rstrip()
        seqr = ref_dict[f]

        fa_handle = open(path.join(args.outdir, f + '.fa'), 'w')
        SeqIO.write(seqr, fa_handle, 'fasta')

        res_name = path.join(args.outdir, f + '.out')
        res_handle = open(res_name, 'w')

        p = subprocess.Popen([ args.iupred, res_name, "long" ], stdout=res_handle)
        p.wait()

        res_handle.close()
        res_handle = open(res_name)

        out_handle = open(path.join(args.outdir, f + '.tab'), 'w')
        for l in res_handle:
            if l.startswith('#'):
                continue

            f = l.strip().split(' ')
            print >>out_handle, '\t'.join([f[0], f[6]])
