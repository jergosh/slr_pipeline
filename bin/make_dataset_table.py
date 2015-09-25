import glob
import argparse
from os import path

import utils


argparser = argparse.ArgumentParser()

argparser.add_argument('--indir', metavar="dir", type=str, required=True)
argparser.add_argument('--outfile', metavar="output_file", type=str, required=True)

def main():
    args = argparser.parse_args()
    outfile = open(args.outfile, 'w')

    for infile in glob.glob(path.join(args.indir, '*', '*_*_*_matched.res')):
        dir, basename = path.split(infile)
        fields = basename.split('_')
        print >>outfile, '\t'.join([ fields[0], '_'.join(fields[1:3]) ])

if __name__ == "__main__":
    main()
