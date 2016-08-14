import os 
from os import path as path
import shutil
import glob
import re
import utils
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import ete2

def match_ids(tree, aln):
    tree_ids = set()

    for l in tree.get_leaves():
        new_label = l.name.partition(' ')[0].upper()
        tree_ids.add(new_label)
        l.name = new_label
    
    aln_ids = set()
    for sr in aln:
        new_id = sr.id.partition(' ')[0].upper()

        aln_ids.add(new_id)

    matched_ids = tree_ids.intersection(aln_ids)

    return matched_ids


fubar_str = """inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]="New Analysis";
inputRedirect["03"]="%s";
inputRedirect["04"]="Default";
inputRedirect["05"]="%s";
inputRedirect["06"]="%s";
inputRedirect["07"]="Neutral";
inputRedirect["08"]="MEME";
inputRedirect["09"]="-1";
inputRedirect["10"]="0.05";
inputRedirect["11"]="N";
inputRedirect["12"]="%s";

ExecuteAFile ("/nfs/research2/goldman/gregs/sw/lib/hyphy/TemplateBatchFiles/QuickSelectionDetection.bf", inputRedirect);
"""

argparser = argparse.ArgumentParser()
argparser.add_argument('--alndir', metavar='input_directory', type=str, required=True)
argparser.add_argument('--treedir', metavar='tree_directory', type=str, required=True)
argparser.add_argument('--filter', metavar='glob', type=str, default="*")
argparser.add_argument('--outdir', metavar='output_directory', type=str, required=True)
argparser.add_argument('--clade', metavar='clade', type=str, required=True)
args = argparser.parse_args()

utils.check_dir(args.outdir)
utils.check_dir(path.join(args.outdir, args.clade))

for input_file in glob.glob(path.join(args.alndir, args.clade, "*", args.filter)):
    print input_file
    basename = path.basename(input_file)
    prefix = basename.partition('_')[0][:2]
    input_name = basename.rpartition('_')[0]
    tree_file = path.join(args.treedir, args.clade, prefix, input_name+'.nh')

    aln = list(SeqIO.parse(open(input_file), 'fasta'))
    tree = ete2.Tree(tree_file)
    matched_ids = match_ids(tree, aln)
    tree.prune(matched_ids)
    filtered_aln = [ SeqRecord(sr.seq, id=sr.id, description="")
                       for sr in aln if sr.id in matched_ids ]

    utils.check_dir(path.join(args.outdir, args.clade, prefix))
    out_fasta = path.abspath(path.join(args.outdir, args.clade, prefix, input_name+'.fa'))
    out_tree = path.abspath(path.join(args.outdir, args.clade, prefix, input_name+'.nwk'))
    ntfit = path.abspath(path.join(args.outdir, args.clade, prefix, input_name+'.ntfit'))

    SeqIO.write(filtered_aln, out_fasta, 'fasta')
    tree.write(format=9, outfile=out_tree)
    control_file = open(path.join(args.outdir, args.clade, prefix, input_name + '_HYPHY'), 'w')
    out_file = path.abspath(path.join(args.outdir, args.clade, prefix, input_name+'.txt'))
    print >>control_file, fubar_str % (out_fasta, out_tree, ntfit, out_file)
