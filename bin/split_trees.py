import os
from os import path
import sys
import pickle
import operator
import time
import re
from pprint import pprint
from argparse import ArgumentParser

import ete2
import emf
import utils

ens_RE = re.compile("[A-Z]+")

argparser = ArgumentParser()

argparser.add_argument('--emf', metavar='emf_file', type=str, required=True)
argparser.add_argument('--clade', metavar='clade_name', type=str, required=True)
argparser.add_argument('--treeroot', metavar='tree_root', type=str, required=True)
argparser.add_argument('--imgroot', metavar='img_root', type=str, required=True)
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
argparser.add_argument('--species_cache', metavar='cache_file', type=str, required=True)
argparser.add_argument('--thr', metavar='value', type=float, default=0.6)


def filter_clades(species, clade_names):
    clade_names = set(clade_names)
    clades = { cl: [] for cl in clade_names }

    parent = None
    for sp in species:
        sp_clades = clade_names.intersection(get_clades(sp))
        for cl in sp_clades:
            clades[cl].append(sp)

    return clades

clade_cache = {}
def get_clades(sp):
    parent = sp
    clades = []
    while parent != "33154": # Eukaryota
        if parent not in clade_cache:
            if parent == "canis familiaris":
                parent = "canis lupus familiaris"
            sp_info = utils.ens_get("taxonomy/id/", parent)
            clade_cache[parent] = sp_info
            time.sleep(0.33)
        else:
            print parent, "cached!"
            sp_info = clade_cache[parent]
        # pprint(sp_info["parent"])
        clades.append(sp_info["name"])
        parent = sp_info["parent"]["id"]

    return clades

def load_clade(fh):
    for l in fh:
        name, taxid = l.rstrip().split('\t')
        self.thr = thr

class TCC(object):
    def __init__(self, clade, thr):
        self.clade = clade
        # self.condition = condition
        self.thr = thr

    def calc(self, node):
        desc_species = set([ n.species for n in node.get_leaves() ])
        clade_cov = float(len(desc_species.intersection(self.clade))) / len(self.clade)

        return clade_cov

    def check(self, node):
        return self.calc(node) > self.thr

class TCCList(object):
    def __init__(self):
        self.tccs = []

    def add(self, tcc):
        self.tccs.append(tcc)

    def check(self, node):
        for tcc in self.tccs:
            desc_species = set([ n.species for n in node.get_leaves() ])
            clade_cov = float(len(desc_species.intersection(tcc.clade))) / len(tcc.clade)
            if not tcc.condition(clade_cov, tcc.thr):
                return False

        return True

def split_tree(intree, tcc):
    print "-"*40
    seqsets, subtrees = [], []

    trees = []
    to_remove = []
    for node in intree.traverse("postorder"):
        if node.dist > 100:
            node.detach()
            print "Detaching:", str(node), len(node.get_leaves())
            if len(node.get_leaves()) > 1:
                trees.append(node)
        # TODO: Does this do anything?
        elif node.is_leaf() and not getattr(node, 'species', None):
            print >>sys.stderr, "REMOVING b/c of lack of species ann.", node.name
            node.detach()
        elif node.is_leaf() and not node.name.startswith("ENS"):
            # print >>sys.stderr, "REMOVING apparently non-Ensembl species ID", node.name
            node.detach()

    if len(intree.get_children()):
        trees.append(intree)

    # if len(to_remove):
    #     print to_remove
    #     for n in to_remove:
    #         n.detach()

    # TMP
    for node in intree.get_leaves():
        if not getattr(node, 'species', None):
            print "MISSING SPECIES", node.name
            print node.get_children()

    trees_fixed = []
    for tree in trees:
        root_to_fix = False
        for node in tree.traverse("postorder"):
            children = node.get_children()
            if len(children) == 1:
                parent = node.up
                if parent is None: # Root has one child
                    print >>sys.stderr, "Fixing root"
                    root_to_fix = True
                    break

                child = node.children[0]

                # print "Reattaching", child.name
                brlen = child.dist + node.dist
                node.detach()
                parent.add_child(child, dist=brlen)
                    
                # print len(parent.get_children())
                # for child in parent.get_children():
                #     print "\t"+str(len(child.get_children()))

        if root_to_fix:
            assert len(tree.get_children()) == 1
            trees_fixed.append(tree.get_children()[0])
        else:
            trees_fixed.append(tree)

    for tree in trees_fixed:
        for node in tree.traverse("postorder"):
            children = node.get_children()

            if not len(children):
                if node.name.startswith("ENSP0"):
                    node.add_features(n_human=1)
                else:
                    node.add_features(n_human=0)

                node.add_features(done=False, split=False)
                continue

            else:
                node.add_features(n_human=(children[0].n_human+children[1].n_human))

            species_codes = []
            for n in node.get_leaves():
                match = ens_RE.match(n.name)
                species_codes.append(match.group())

            paralog_frac = float(len(species_codes)-len(set(species_codes)))/len(species_codes)
            node.add_features(paralog_frac=paralog_frac)
            node.add_features(tcc=tcc.calc(node))
            node.add_features(split=tcc.check(node))

            if any([ n.done for n in children ]):
                node.add_features(done=True)
                for child in children:
                    if child.split and not child.done:
                        seqsets.append([ (n.name, n.species) for n in child.get_leaves() ])
                        subtrees.append(child)

                continue

            if len(children) != 2:
                print node
                print children
                sys.exit(-1)

            # if any([ child.split for child in children ]):
            #     node.done = True
            # else:
            #     node.done = False

            # if all([ child.split for child in children ]):
            #     node.done = True
            # else:
            #     node.done =False
                
            if (children[0].split and children[1].split) or node.is_root():
                node.add_features(done=True)

            elif children[0].split or children[1].split:
                if children[0].split:
                    child_split = children[0]
                    child_unsplit = children[1]
                else:
                    child_split = children[1]
                    child_unsplit = children[0]

                if child_split.n_human and child_unsplit.n_human:
                    if node.paralog_frac <= child_split.paralog_frac:
                        node.add_features(done=False)
                    else:
                        node.add_features(done=True)
                elif child_split.n_human and not child_unsplit.n_human:
                    if node.paralog_frac <= child_split.paralog_frac:
                        node.add_features(done=False)
                    else:
                        node.add_features(done=True)
                elif not child_split.n_human and child_unsplit.n_human:
                    # Continuing up the tree is the only chance for this human ID
                    # to be included
                    node.add_features(done=False)
                # Neither subtree contains a human sequence -- we continue up the tree
                elif not child_split.n_human and not child_unsplit.n_human:
                    node.add_features(done=False)
                    # Is it 'fair' to include 
                
            else:
                node.add_features(done=False)

            if node.done:
                for child in children:
                    if child.split:
                        seqsets.append([ (n.name, n.species) for n in child.get_leaves() ])
                        subtrees.append(child)

    return seqsets, subtrees

colours = { -1: "black", 0: "red", 1: "green", 2: "blue", 3: "purple", 4: "orange", 5: "yellow", 6: "grey", 7: "#009999",
            8: "#FF7400"}
def make_layout(nodesets):
    nodemap = {}
    for i, ns in enumerate(nodesets):
        for n in ns:
            nodemap[n[0]] = i

    def layout(node):
        if node.D == "Y":
            node.img_style["shape"] = "circle"
            node.img_style["size"] = 12
            node.img_style["fgcolor"] = "red"

        if node.is_leaf():
            nameFace = ete2.faces.AttrFace("name", fsize=20, 
                                           fgcolor=colours[nodemap.get(node.name, -1)])
            ete2.faces.add_face_to_node(nameFace, node, column=0)

    return layout

def main():
    args = argparser.parse_args()

    emf_file = args.emf
    out_root = args.outroot
    img_root = args.imgroot
    tree_root = args.treeroot
    clades_pickle = args.species_cache

    all_species = utils.ens_get("/info/species/")["species"]
    all_species_names = [ it["name"].replace("_", " ") for it in all_species ]
    # FIXME Temporary
    # all_species_names.remove("Ancestral sequences")

    if path.exists(clades_pickle):
        Clades = pickle.load(open(clades_pickle, 'rb'))
    else:
        Clades = filter_clades(all_species_names,
                               [ "Eutheria", "Glires", "Laurasiatheria", "Sauria", "Mammalia", "Primates" ])
        pickle.dump(Clades, open(clades_pickle, 'wb'))

    pprint(Clades)

    # TL = TCCList()
    # TL.add(TCC(Clades[args.clade], operator.ge, args.thr))
    tcc = TCC(Clades[args.clade], args.thr)

    utils.check_dir(path.join(out_root, args.clade))

    tree_id = 1
    for tree in emf.EMF(emf_file):
        print tree_id
        treedir = path.join(tree_root, str(tree_id)[:2])
        utils.check_dir(treedir)

        tree.write(outfile=path.join(treedir, "{}.nh".format(tree_id)))

        seqsets, subtrees = split_tree(tree, tcc)
        outdir = path.join(out_root, args.clade, str(tree_id)[:2])
        utils.check_dir(outdir)

        # Treevis
        # layout = make_layout(seqsets)
        # imgdir = path.join(img_root, args.clade)
        # utils.check_dir(imgdir)
        # imgfile = path.join(imgdir, "{}.pdf".format(tree_id))
        # tree.render(imgfile, layout=layout)

        set_id = 1
        for seqset, subtree in zip(seqsets, subtrees):
            outfile = open(path.join(outdir, "{}_{}.tab".format(tree_id, set_id)), 'w')
            for seqid in seqset:
                print >>outfile, '\t'.join(seqid)

            subtree.write(outfile=path.join(outdir, "{}_{}.nh".format(tree_id, set_id)), 
                          format=6)
            set_id += 1

        tree_id += 1

if __name__ == "__main__":
    main()
