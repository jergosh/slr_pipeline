from glob import glob
import os
from os import path
import argparse
import re

from dendropy import Tree
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def match_ids(tree, fasta):
    tree_ids = set([ l.taxon.label for l in tree.leaf_nodes() ])
    fasta_ids = set([ sr.id for sr in fasta ])
    matched_ids = tree_ids.intersection(fasta_ids)

    return matched_ids

def write_paml(aln, fh):
    print >>fh, len(aln), aln.get_alignment_length()
    for seqr in aln:
        print >>fh, seqr.id
        print >>fh, seqr.seq

def write_slr(slr_fh, nt_aln_fh, tree_file, gene_name):
    slr_template = """seqfile: %s
treefile: %s
outfile: %s.res
positive_only: 0"""
    print >>slr_fh, slr_template % (path.basename(nt_aln_fh.name), path.basename(tree_file.name),
                                    path.basename(gene_name))



argparser = argparse.ArgumentParser()

argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--treeroot', metavar='tree_root', type=str, required=True)
argparser.add_argument('--alnroot', metavar='aln_root', type=str, required=True)
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
argparser.add_argument('--logroot', metavar='log_root', type=str, required=True)

args = argparser.parse_args()

treeroot = args.treeroot
alndir = args.alnroot
slrdir = args.outroot
logroot = args.logroot

utils.check_dir(slrdir)
utils.check_dir(logroot)

# for infile in glob(path.join(alndir, args.clade, "*", "*_prank.best.fas")):
for infile in glob(path.join(alndir, args.clade, "*", "*_prank.best.fas"))[:2]:
    basename = path.basename(infile).rpartition('_')[0]

    treedir = path.join(treeroot, args.clade, basename[:2])
    treefile = path.join(treedir, basename + '.nh')

    outdir = path.join(slrdir, args.clade, basename[:2])
    utils.check_dir(outdir)

    fasta = [ f for f in SeqIO.parse(open(infile), 'fasta') ]
    tree = Tree.get_from_path(treefile, 'newick')

    matched_ids = match_ids(tree, fasta)
    tree.retain_taxa_with_labels(matched_ids)
    tree.deroot()

    aln_ids = {}
    for idx, seqr in enumerate(fasta):
        # Sort out the ID for the tree
        aln_ids[seqr.id] = str(idx+1)

    for node in tree:
        if node.edge.length == 0:
            node.edge.length = 0.00001
        if node.is_leaf():
            node.label = aln_ids[node.taxon.label]
        else:
            node.label = ""


    pamlfile = open(path.join(outdir, basename + '_slr.paml'), 'w')
    outtree = open(path.join(outdir, basename + '_slr.nwk'), 'w')

    print >>outtree, len(fasta), "1"
    print >>outtree, tree.as_string('newick', suppress_rooting=True, 
                                    suppress_internal_taxon_labels=True,
                                    suppres_internal_node_labels=True)

    write_paml(MultipleSeqAlignment(fasta), pamlfile)

    ctl_file = file(path.join(outdir, basename + '_matched.ctl'), 'w')
    write_slr(ctl_file, pamlfile, outtree, path.join(outdir, basename + '_matched'))
        
