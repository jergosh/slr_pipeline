import glob
import argparse
import sys
from os import path

import ete2
from Bio import SeqIO

import utils

argparser = argparse.ArgumentParser()

# Use an ID list + dataset map as a way of limiting the number of files?
argparser.add_argument('--treedir', metavar="dir", type=str, required=True)
# argparser.add_argument('--alndir', metavar="dir", type=str, required=True)
argparser.add_argument('--domaindir', metavar="dir", type=str, required=True)
argparser.add_argument('--dataset_map', metavar="file", type=str, required=True)
argparser.add_argument('--cath_map', metavar="file", type=str, required=True)
argparser.add_argument('--outdir', metavar="output_file", type=str, required=True)

raxml_cmd = "raxmlHPC -T 2 -t {} -s {} -n {} -f e -p 1 -m PROTEINGAMMA" # tree, aln, name of the run

# The trees are in seqsets
# We might need to prune to the sequences present in the alignment
# Run RAxML on each of them in the queue

def main():
    args = argparser.parse_args()

    utils.check_dir(args.outdir)

    dataset_map = {}
    for l in open(args.dataset_map):
        f = l.rstrip().split('\t')
        dataset_map[f[0]] = f[1]
  
    for l in open(args.cath_map):
        stable_id, ranges, cath_id, pdb_id = l.rstrip().split('\t')
        dataset = dataset_map[stable_id]
        prefix = dataset.partition('_')[0][:2]
        intree = ete2.Tree(path.join(args.treedir, prefix, dataset+'.nh'))
        aln_fn = path.absolute(path.join(args.domaindir, prefix, dataset+'.fa'))
        aln = SeqIO.to_dict(SeqIO.parse(aln_fn, 'fasta'))
        outtree_fn = path.join(args.outdir, prefix, dataset+'.nh')
        
        utils.check_dir(path.join(args.outdir, prefix))

        intree.prune(aln.keys())
        intree.write(outfile=outtree_fn, format=9)

        logfile = path.join(args.outdir, prefix, dataset+'.log')

        raxml = raxml_cmd.format(tree_fn, aln_fn, dataset)
        p = Popen(['bsub', '-o'+logfile, '-n 2', # '-g', '/slr',
                   '-cwd', args.outdir, raxml])
        p.wait()


if __name__ == "__main__":
    main()