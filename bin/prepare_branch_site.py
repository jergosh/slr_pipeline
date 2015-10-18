import argparse
from os import path
import glob
# from Bio.Phylo.PAML import codeml

# Model A (alternative): model = 2, NSsites = 2,  fix_omega = 0 
# Model A1 (null): model = 2, NSsites = 2,  fix_omega = 1, omega = 1 
from Bio import AlignIO
from Bio import SeqIO
# from Bio.PAML import codeml
import ete2

import utils

argparser = argparse.ArgumentParser()

argparser.add_argument('--infile', metavar="dir", type=str, required=True)
argparser.add_argument('--dataset_map', metavar="dir", type=str, required=True)
argparser.add_argument('--treeroot', metavar="dir", type=str, required=True)
argparser.add_argument('--alnroot', metavar="dir", type=str, required=True)
argparser.add_argument('--outroot', metavar="dir", type=str, required=True)

# Alignment stays the same but maybe needs to be copied to a subdir
# Dir structure:
# 100 dirs 
# subdir for each dataset/branch combination, in the form 100_1_2
def main():
    args = argparser.parse_args()

    dataset_map = {}
    for l in open(args.dataset_map):
        f = l.rstrip().split('\t')
        dataset_map[f[0]] = f[1]

    utils.check_dir(args.outroot)

    processed = set()
    for l in open(args.infile):
        f = l.rstrip().split('\t')
        stable_id = f[0]
        if stable_id in processed:
            continue
        else:
            processed.add(stable_id)

        dataset = dataset_map[stable_id]
        basename = dataset+'.nh'
        prefix = dataset.partition('_')[0][:2]
        tree_fn = path.join(args.treeroot, prefix, basename)
        print "Processing", dataset

        tree = ete2.Tree(tree_fn)
        tree.unroot()
        aln = AlignIO.read(path.join(args.alnroot, prefix, dataset+"_prank.best.fas"), 'fasta')
        seqnames = SeqIO.to_dict(aln).keys()

        outdir = path.join(args.outroot, prefix)
        dataset_dir = path.join(outdir, dataset)
        utils.check_dir(outdir)
        utils.check_dir(dataset_dir)

        tree.prune(seqnames)
        utils.write_paml(open(path.join(outdir, dataset+'.paml'), 'w'), aln)

        for i, node in enumerate(tree.traverse("postorder")):
            old_name = node.name
            if node.name == 'NoName':
                node.name = '#1'
            else:
                node.name += ' #1'

            subset_dir = path.join(dataset_dir, dataset+'_'+str(i))
            utils.check_dir(subset_dir)
            tree.write(outfile=path.join(subset_dir, basename), format=3)

            for run_id in [ "1", "2" ]:
                workingdir = path.join(subset_dir, run_id)
                utils.check_dir(workingdir)
                # Probably the null can be run just once for all the branches
            node.name = old_name
        

if __name__ == "__main__":
    main()
