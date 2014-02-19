from dendropy import Tree
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser

import re
sp_RE = re.compile("([A-Z]+)")

argparser = ArgumentParser()

argparser.add_argument('--fasta', metavar='fasta_file', type=str, required=True)
argparser.add_argument('--tree', metavar='tree_file', type=str, required=True)
argparser.add_argument('--outroot', metavar='outfile_root', type=str, required=True)

def match_ids(tree, fasta):
    tree_ids = set([ l.taxon.label for l in tree.leaf_nodes() ])
    fasta_ids = set([ sr.id for sr in fasta ])
    matched_ids = tree_ids.intersection(fasta_ids)

    return matched_ids

if __name__ == "__main__":
    args = argparser.parse_args()

    output_fasta = args.outroot + ".fa"
    output_tree = args.outroot + ".nh"

    fasta = [ f for f in SeqIO.parse(open(args.fasta), 'fasta') ]
    tree = Tree.get_from_path(args.tree, 'newick')

    matched_ids = match_ids(tree, fasta)
    tree.retain_taxa_with_labels(matched_ids)

    for node in tree:
        if node.is_leaf():
            # pass
            node.taxon.label = node.taxon.label
        else:
            node.label = ""

    # filtered_fasta = [ sr for sr in fasta if sr.id in matched_ids ]
    filtered_fasta = [ SeqRecord(sr.seq, id=sr.id, description="")
                       for sr in fasta if sr.id in matched_ids ]

    SeqIO.write(filtered_fasta, output_fasta, 'fasta')
    tree.write_to_path(output_tree, 'newick', suppress_rooting=True)
