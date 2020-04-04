import sys
import glob
import pickle
import argparse 
from os import path

import PAML

argparser = argparse.ArgumentParser()

argparser.add_argument('--inroot', metavar="dir", type=str, required=True)
argparser.add_argument('--dataset_map', metavar="dir", type=str, required=True)
# argparser.add_argument('--treeroot', metavar="dir", type=str, required=True)
# argparser.add_argument('--alnroot', metavar="dir", type=str, required=True)
# argparser.add_argument('--outroot', metavar="dir", type=str, required=True)

# TODO Output two or three files 
# Per-gene-branch combination -- with the LRT
# Per gene-branch-site -- with the categories/posteriors
# Only with positives sites -- equivalently we could just count the sites

def main():
    args = argparser.parse_args()

    dataset_map = {}
    for l in open(args.dataset_map):
        f = l.rstrip().split('\t')
        dataset_map[f[0]] = f[1]

    for indir in glob.glob(path.join(args.inroot, '*', '*', '*[0-9]')):
        basename = path.basename(indir)
        f = basename.split('_')
        dataset = '_'.join(f[:2])
        sample = f[2]

        pickle_1 = path.join(indir, '1', "results_1.pk")
        pickle_2 = path.join(indir, '2', "results_2.pk")
        if not path.exists(pickle_1) or not path.exists(pickle_2):
            # print >>sys.stderr, "Skipping", basename
            continue

        lnL_1 = pickle.load(open(pickle_1))["NSsites"][2]["lnL"]
        lnL_2 = pickle.load(open(pickle_2))["NSsites"][2]["lnL"]
        # print dataset, sample, 2*(lnL_2 - lnL_1)

        rstfile = path.join(indir, "2", "rst")
        for site in PAML.rstparser(open(rstfile)):
            if site[5] > 0.5:
                print '\t'.join([ str(i) for i in [ dataset, sample, 2*(lnL_2 - lnL_1) ] + site])

if __name__ == "__main__":
    main()
