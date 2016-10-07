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
argparser.add_argument('--iupred', metavar="iupred", type=str, default="iupred")


if __name__ == "__main__":
    args = argparser.parse_args()

    ref_handle = open(args.reffile)
    ref_dict = SeqIO.to_dict(SeqIO.parse(ref_handle, "fasta"))

    infile = open(args.infile)

    for l in infile:
        ens_id = l.rstrip()
        seqr = ref_dict[f]
        seqr.description = ""
        
        fa_handle = open(path.join(args.outdir, ens_id + '.fa'), 'w')
        SeqIO.write(seqr, fa_handle, 'fasta')
        fa_handle.close()
        
        res_name = path.join(args.outdir, ens_id + '.out')
        res_handle = open(res_name, 'w')

        p = subprocess.Popen([ args.iupred, fa_handle.name, "long" ], stdout=res_handle,
                             env={ "IUPred_PATH": path.dirname(args.iupred) })
        p.wait()

        res_handle.close()
        res_handle = open(res_name)

        out_handle = open(path.join(args.outdir, ens_id + '.tab'), 'w')
        for r in res_handle:
            if r.startswith('#'):
                continue

            f = r.strip().split(' ')
            print >>out_handle, '\t'.join([f, f[0], f[6]])
