import sys
from os import path
import glob
import re
import argparse
from subprocess import Popen

from Bio import SeqIO
import utils

ens_RE = re.compile("[A-Z]*")
cmd_template = "/nfs/research2/goldman/gregs/sw/SCRATCH-1D_1.0/bin/run_SCRATCH-1D_predictors.sh {} {}"

argparser = argparse.ArgumentParser()

argparser.add_argument('--indir', metavar="input_dir", type=str, required=True)
argparser.add_argument('--outdir', metavar="output_dir", type=str, required=True)
argparser.add_argument('--refgenome', metavar="output_dir", type=str, required=True)

def main():
    args = argparser.parse_args()

    ref = SeqIO.to_dict(SeqIO.parse(open(args.refgenome), 'fasta'))

    utils.check_dir(args.outdir)
    for f in glob.glob(path.join(args.indir, "*", "*.fa")):
        for seqr in SeqIO.parse(open(f), 'fasta'):
            species_code = ens_RE.match(seqr.id).group(0)
            if species_code == "ENSP":
                seqr = ref[seqr.id]

                indir, basename = path.split(f)
                outroot, subdir = path.split(indir)
                input_file = path.join(args.outdir, subdir, seqr.id+'_'+basename) 
                log_file = path.join(args.outdir, subdir, seqr.id+'_'+basename+'.log')

                print input_file
                utils.check_dir(path.join(args.outdir, subdir))
                SeqIO.write([ seqr ], open(input_file, 'w'), 'fasta')
                cmd = cmd_template.format(input_file, path.join(args.outdir, subdir, seqr.id+'_'+basename))
                p = Popen(['bsub', '-o'+log_file, cmd])
        
                p.wait()


if __name__ == "__main__":
    main()
