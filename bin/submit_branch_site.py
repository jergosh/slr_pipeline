import argparse
from os import path
import glob
import utils
from subprocess import Popen

import ete2

bsub_cmd = "python bin/paml_bs_wrapper.py --indir {} --dataset {} --sample {}"

argparser = argparse.ArgumentParser()

argparser.add_argument('--infile', metavar="dir", type=str, required=True)
argparser.add_argument('--dataset_map', metavar="dir", type=str, required=True)
argparser.add_argument('--indir', metavar="dir", type=str, required=True)
argparser.add_argument('--rerun', metavar="bool", type=bool, default=False)

def main():
    args = argparser.parse_args()

    dataset_map = {}
    for l in open(args.dataset_map):
        f = l.rstrip().split('\t')
        dataset_map[f[0]] = f[1]

    processed = set()
    for l in open(args.infile):
        f = l.rstrip().split('\t')
        stable_id = f[0]
        if stable_id in processed:
            continue
        else:
            processed.add(stable_id)

        dataset = dataset_map[stable_id]

        prefix = dataset.partition('_')[0][:2]
        outdir = path.join(args.indir, prefix)
        dataset_dir = path.join(outdir, dataset)
        basename = dataset+'.nh'
        # tree_fn = path.join(dataset_dir, basename)

        # tree = ete2.Tree(tree_fn)

        for sample_dir in glob.glob(path.join(dataset_dir, "*[0-9]")):
            i = path.basename(sample_dir).rpartition('_')[2]
            subset_dir = path.join(dataset_dir, dataset+'_'+i)
            
            if args.rerun and path.exists(path.join(subset_dir, "1", "results_1.pk")) and \
                    path.exists(path.join(subset_dir, "2", "results_2.pk")):
                print "Skipping", subset_dir
                continue

            logfile = path.join(subset_dir, 'paml.log')

            bsub = bsub_cmd.format(subset_dir, dataset, i)
            p = Popen(['bsub', '-o'+logfile, bsub])
            p.wait()


if __name__ == "__main__":
    main()
