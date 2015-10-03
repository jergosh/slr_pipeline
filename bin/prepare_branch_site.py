import argparse
from os import path
import glob
# from Bio.Phylo.PAML import codeml

# Model A (alternative): model = 2, NSsites = 2,  fix_omega = 0 
# Model A1 (null): model = 2, NSsites = 2,  fix_omega = 1, omega = 1 

import ete2

argparser = argparse.ArgumentParser()

# Use an ID list + dataset map as a way of limiting the number of files?
argparser.add_argument('--treeroot', metavar="dir", type=str, required=True)

# Alignment stays the same but maybe needs to be copied to a subdir
# Dir structure:
# 100 dirs 
# subdir for each dataset/branch combination, in the form 100_1_2
def main():
    args = argparser.parse_args()

    for tree_fn in glob.glob(path.join(args.treeroot, '*', '*.nh')):
        basename = path.basename(tree_fn)
        dataset = basename.partition('_')[0]
        prefix = dataset.partition('_')[0][:2]

        tree = ete2.Tree(tree_fn)

        for node in tree.traverse("postorder"):
            print node.name
            node.name += "#1"

if __name__ == "__main__":
    main()
