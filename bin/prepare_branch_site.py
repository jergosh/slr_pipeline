import argparse
from os import path
import glob
# from Bio.Phylo.PAML import codeml

# Model A (alternative): model = 2, NSsites = 2,  fix_omega = 0 
# Model A1 (null): model = 2, NSsites = 2,  fix_omega = 1, omega = 1 

import ete2

import utils

argparser = argparse.ArgumentParser()

# Use an ID list + dataset map as a way of limiting the number of files?
argparser.add_argument('--treeroot', metavar="dir", type=str, required=True)
argparser.add_argument('--outroot', metavar="dir", type=str, required=True)

# Alignment stays the same but maybe needs to be copied to a subdir
# Dir structure:
# 100 dirs 
# subdir for each dataset/branch combination, in the form 100_1_2
def main():
    args = argparser.parse_args()

    utils.check_dir(args.outroot)

    for tree_fn in glob.glob(path.join(args.treeroot, '*', '*.nh')):
        basename = path.basename(tree_fn)
        dataset = basename.partition('.')[0]
        prefix = dataset.partition('_')[0][:2]

        tree = ete2.Tree(tree_fn)
        outdir = path.join(args.outroot, prefix)
        utils.check_dir(outdir)
        for i, node in enumerate(tree.traverse("postorder")):
            old_name = node.name
            if node.name == 'NoName':
                node.name = '#1'
            else:
                node.name += ' #1'

            workingdir = path.join(outdir, dataset+'_'+str(i))
            utils.check_dir(workingdir)
            tree.write(outfile=path.join(workingdir, basename))
            node.name = old_name
            print workingdir, node.name
        

if __name__ == "__main__":
    main()
