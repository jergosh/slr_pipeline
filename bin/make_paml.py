from os import path
from dendropy import Tree
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from argparse import ArgumentParser

from utils import write_paml

def write_slr(slr_fh, nt_aln_fh, tree_file, gene_name):
    slr_template = """seqfile: %s
treefile: %s
outfile: %s.res
positive_only: 0"""
    print >>slr_fh, slr_template % (path.basename(nt_aln_fh.name), path.basename(tree_file.name),
                                    path.basename(gene_name))

argparser = ArgumentParser()

argparser.add_argument('--fasta', metavar='fasta_file', type=str, required=True)
argparser.add_argument('--tree', metavar='tree_file', type=str, required=True)
argparser.add_argument('--outroot', metavar='outfile_root', type=str, required=True)
argparser.add_argument('--paml',dest='paml',action='store_true')
argparser.set_defaults(paml=False)

args = argparser.parse_args()

tree = Tree.get_from_path(args.tree, 'newick')
aln = AlignIO.read(open(args.fasta), 'fasta')

if args.paml:
    paml_file = open(args.outroot + '.paml', 'w')
    tree_file = open(args.outroot + '.nwk', 'w')
else:
    paml_file = open(args.outroot + '_slr.paml', 'w')
    tree_file = open(args.outroot + '_slr.nwk', 'w')


aln_ids = {}
for idx, seqr in enumerate(aln):
    # Sort out the ID for the tree
    aln_ids[seqr.id] = str(idx+1)


for node in tree:
    if not args.paml:
        if node.is_leaf():
            node.label = aln_ids[node.taxon.label]
        else:
            node.label = ""

print >>tree_file, len(aln), "1"
print >>tree_file, tree.as_string('newick', suppress_rooting=True, 
                                  suppress_edge_lengths=args.paml,
                                  suppress_internal_taxon_labels=True,
                                  suppres_internal_node_labels=True)

write_paml(MultipleSeqAlignment(aln), paml_file)

ctl_file = file(args.outroot + '.ctl', 'w')
write_slr(ctl_file, paml_file, tree_file, args.outroot)
