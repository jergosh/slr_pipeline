import os
from os import path
import pandas
import sys
import argparse
from subprocess import Popen
import glob 

preamble = """
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
  <head>
    <meta http-equiv="content-type" content="text/html; charset=utf-8"/>
    <title>title</title>
    <link rel="stylesheet" type="text/css" href="style.css"/>
    <script type="text/javascript" src="script.js"></script>
  </head>
  <body>
"""

postamble = """
</body>
</html>
"""


def main():
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--resdir', metavar='result_dir', type=str, required=True)
    argparser.add_argument('--treedir', metavar='tree_dir', type=str, required=True)
    argparser.add_argument('--alndir', metavar='aln_dir', type=str, required=True)
    argparser.add_argument('--signature', metavar='filename', type=str, default="ENS*.res")
    argparser.add_argument('--outalns', metavar='output_aln', type=str, required=True)
    argparser.add_argument('--outtrees', metavar='output_tree', type=str, required=True)

    args = argparser.parse_args()

    tree_content = [ preamble ]
    aln_content = [ preamble ]
    for resfile in glob.glob(path.join(args.resdir, '*', args.signature)):
        basename = path.basename(resfile).rpartition('.')[0]
        fields = basename.split('_')
        ens = fields[0]
        dataset = '_'.join(fields[1:3])
        print ens, dataset

        tree = path.join("trees", fields[1][:2], dataset+'.nh')
        aln = path.join("aln", fields[1][:2], dataset+'_prank.best.fas') 
        tree_link = "<a href=\"" + tree + "\">" + ens + "</a>"
        aln_link = "<a href=\"" + aln + "\">" + ens + "</a>"

        tree_content.append("<p>"+tree_link+"</p>")
        aln_content.append("<p>"+aln_link+"</p>")

    tree_content.append(postamble)
    aln_content.append(postamble)

    tree_content = '\n'.join(tree_content)
    aln_content = '\n'.join(aln_content)

    tree_file = open(args.outtrees, 'w')
    aln_file = open(args.outalns, 'w')

    print >>tree_file, tree_content
    print >>aln_file, aln_content

if __name__ == "__main__":
    main()
