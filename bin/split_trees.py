import os
from os import path
import sys
import pickle
import operator
import time
import re
import copy
from pprint import pprint
from argparse import ArgumentParser

import ete3
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
argparser.add_argument('--species_list', metavar='species_file', type=str, required=True)
argparser.add_argument('--thr', metavar='value', type=float, default=0.6)

def load_species(fh):
    species_names, species = zip(*[ l.rstrip().split('\t')[1:] for l in fh ])
    return [ s.lower() for s in species_names ], species

def remove_node(node):
    parent = node.up
    if parent is None:
        return node 

    grandparent = parent.up 

    node.detach()
    assert len(parent.children) == 1
    sibling = parent.children[0].detach()
    dist = parent.dist
    parent.detach()

    if grandparent is not None:
        grandparent.add_child(sibling)
        assert sibling.up == grandparent
        sibling.dist += dist

    return node

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
        if len(node.children) == 0:
            return 0.0

        desc_species = set([ n.species for n in node.get_leaves() ])
        clade_cov = float(len(desc_species)) / len(self.clade)

        return clade_cov

    def check(self, node):
        return self.calc(node) >= self.thr

class TCC_abs(TCC):
    def __init__(self, clade, thr):
        super(TCC_abs, self).__init__(clade, thr)

    def calc(self, node):
        desc_species = set([ n.species for n in node.get_leaves() ])
        clade_cov = float(len(desc_species.intersection(self.clade))) / len(self.clade)

        return clade_cov

    def check(self, node):
        return self.calc(node) >= self.thr

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

abandoned = 0
def split_tree(intree, tcc):
    global abandoned
    seqsets, subtrees = [], []

    for node in intree.traverse("postorder"):
        species_codes = []
        for n in node.get_leaves():
            match = ens_RE.match(n.name)
            species_codes.append(match.group())

        paralog_frac = float(len(species_codes)-len(set(species_codes)))/len(species_codes)

        node.add_features(paralog_frac=paralog_frac)
        node.add_features(tcc=tcc.calc(node))
        node.add_features(split=tcc.check(node))

    for node in copy.deepcopy(intree).traverse("postorder"):
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
                    print "'Abandoning' at TCC", child_unsplit.tcc, "("+str(child_unsplit.n_human)+")"
                    abandoned += 1
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

            if node.is_root():
                if len(children):
                    if not any([ child.split for child in children ]):
                        if node.split:
                            seqsets.append([ (n.name, n.species) for n in node.get_leaves() ])
                            subtrees.append(node)


    return seqsets, subtrees

## What I've learned here is that setting node as done can fail if both children have been detached
def baseline_tcc(intree, tcc):
    seqsets, subtrees = [], []
    for node in intree.traverse("postorder"):
        children = node.get_children()
        node.add_features(split=tcc.check(node))
        node.add_features(tcc=tcc.calc(node))
        node.add_features(done=False)

        if not len(children):
            if node.name.startswith("ENSP0"):
                node.add_features(n_human=1)
            else:
                node.add_features(n_human=0)
        else:
            node.add_features(n_human=(children[0].n_human+children[1].n_human))

    for node in copy.deepcopy(intree).traverse("postorder"):
        if node.up is None and node.done == False:
            seqsets.append([ (n.name, n.species) for n in node.get_leaves() if getattr(n, 'species', None) ])
            subtrees.append(node.detach())

        elif node.n_human == 1 and node.up.n_human > 1:
            node.up.done = True

            seqsets.append([ (n.name, n.species) for n in node.get_leaves() ])
            subtrees.append(node.detach())


    return seqsets, subtrees

def split_tree_gregj(intree, tcc):
    # for node in intree.traverse("postorder"):
    #     node.add_features(split=tcc.check(node))
    #     node.add_features(tcc=tcc.calc(node))
    #     node.add_features(done=False)

    seqsets, subtrees = [], []

    tree = intree # copy.deepcopy(intree)
    for node in tree.traverse("postorder"):
        children = node.get_children()
        if not len(children) == 2:
            if len(children) == 0:
                continue
            else:
                print "Children", children
                print >>sys.stderr, "Something's up!"
                sys.exit(-1)

        # if any([ ch.done for ch in children ]):
        #     print "Already done"
        #     node.done = True
        #     continue
        
        # if children[0].split and children[1].split:
        if tcc.check(children[0]) and tcc.check(children[1]):
            seqset = [ (n.name, n.species) for n in children[0].get_leaves() ]
            if len(seqset) < 21:
                print seqset
                sys.exit(-1)
            seqsets.append(seqset)
            subtrees.append(children[0].detach())
            seqset = [ (n.name, n.species) for n in children[1].get_leaves() ]
            if len(seqset) < 21:
                print seqset
                sys.exit(-1)
            seqsets.append(seqset)
            subtrees.append(children[1].detach())

            parent = node.up
            remove_node(node)
            # node.done = True
        else:
            if node.is_root():
                print "At root"
                # if node.split:
                if tcc.check(node):
                    print node
                    seqset = [ (n.name, n.species) for n in node.get_leaves() ]
                    print len(seqset)
                    if len(seqset) < 21:
                        sys.exit(-1)
                    seqsets.append(seqset)
                    
                    subtrees.append(node)
                    

    return seqsets, subtrees

def split_tree_old(intree, tcc):
    seqsets, subtrees = [], []

    for node in intree.traverse("postorder"):
        children = node.get_children()
        if any([ n.done for n in children ]):
            # if not all([ n.done for n in children ]):
            #     print "Warning!"
            node.add_features(done=True)
            continue

        if len(children) and any([ n.split for n in children ]):
            seqset = []
            for ch in children:
                if ch.split:
                    seqsets.append([ (n.name, n.species) for n in ch.get_leaves() ])
                    subtrees.append(ch)
            node.add_features(done=True)
        else:
            if tcc.check(node):
                node.add_features(split=True)
            else:
                node.add_features(split=False)

            node.add_features(done=False)

    return seqsets, subtrees
    
def process_tree(intree, fix=True):
    trees = []
    to_keep = []
    for l in intree.iter_leaves():
        if getattr(l, 'species', None):
            to_keep.append(l)

    try:
        intree.prune(to_keep)
    except ete3.coretype.tree.TreeError:
        return []

    for node in intree.traverse("postorder"):
        if fix and node.dist > 100:
            print "Fixing..."
            node.detach()
            # print "Detaching:", str(node), len(node.get_leaves())
            if len(node.get_leaves()) > 1:
                trees.append(node)
        # TODO: Does this do anything?
        elif node.is_leaf() and not getattr(node, 'species', None):
            print >>sys.stderr, "REMOVING b/c of lack of species ann.", node.name
            node.detach()
        elif node.is_leaf() and not node.name.startswith("ENS"):
            print >>sys.stderr, "REMOVING apparently non-Ensembl species ID", node.name
            # node.detach()

    if len(intree.get_children()):
        trees.append(intree)

    trees_fixed = []
    for tree in trees:
        root_to_fix = False
        for node in tree.traverse("postorder"):
            children = node.get_children()
            if len(children) == 1:
                parent = node.up
                if parent is None: # Root has one child
                    # print >>sys.stderr, "Fixing root"
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

            while len(tree.get_children()) == 1:
                tree = tree.get_children()[0]

            if len(tree.get_children()) == 0:
                continue

            assert len(tree.get_children()) == 2
            trees_fixed.append(tree)
        else:
            trees_fixed.append(tree)

    return trees_fixed


colours = { -1: "black", 0: "red", 1: "green", 2: "blue", 3: "purple", 4: "orange", 5: "yellow", 6: "grey", 7: "#009999", 8: "#FF7400", 9: "Wheat", 10: "MistyRose", 11: "Goldenrod", 12: "Brown", 13: "DeepPink", 14: "Gold", 15: "Plum", 16: "Magenta", 17: "DarkSlateBlue", 18: "Olive", 19: "LightCoral", 20: "SteelBlue" }
def make_layout(nodesets):
    nodemap = {}
    for i, ns in enumerate(nodesets):
        for n in ns:
            nodemap[n[0]] = i

    def layout(node):
        node.add_face(ete3.AttrFace("tcc", formatter="%0.2f"), column=0, position="branch-right")

        # if node.D == "Y":
        #     node.img_style["shape"] = "circle"
        #     node.img_style["size"] = 12
        #     node.img_style["fgcolor"] = "red"

        if node.is_leaf():
            nameFace = ete3.faces.AttrFace("name", fsize=20, 
                                           fgcolor=colours[nodemap.get(node.name, -1)])
            ete3.faces.add_face_to_node(nameFace, node, column=0)

    return layout

def main():
    args = argparser.parse_args()

    emf_file = args.emf
    out_root = args.outroot
    out_root_ann = args.outroot+'_ann'
    img_root = args.imgroot
    tree_root = args.treeroot
    ann_root = tree_root+'_ann'
    clades_pickle = args.species_cache

    all_species = utils.ens_get("/info/species/")["species"]
    all_species_names = [ it["name"].replace("_", " ") for it in all_species ]
    # FIXME Temporary
    # all_species_names.remove("Ancestral sequences")

    species_names, species = load_species(open(args.species_list))

    # if path.exists(clades_pickle):
    #     Clades = pickle.load(open(clades_pickle, 'rb'))
    # else:
    #     Clades = filter_clades(all_species_names,
    #                            [ "Eutheria", "Glires", "Laurasiatheria", "Sauria", "Mammalia", "Primates" ])
    #     pickle.dump(Clades, open(clades_pickle, 'wb'))

    # pprint(Clades)

    # clade = set(Clades[args.clade]).intersection(species_names)
    clade = species_names
    print len(clade), clade
    # TL = TCCList()
    # TL.add(TCC(Clades[args.clade], operator.ge, args.thr))
    tcc = TCC(clade, args.thr)

    utils.check_dir(path.join(tree_root))
    utils.check_dir(path.join(out_root))
    utils.check_dir(path.join(out_root, args.clade))

    utils.check_dir(path.join(ann_root))
    utils.check_dir(path.join(out_root_ann))
    utils.check_dir(path.join(out_root_ann, args.clade))

    tree_id = 1
    for tree in emf.EMF(emf_file):
        print tree_id, len(tree.get_leaves()),
        anndir = path.join(ann_root, str(tree_id)[:2])
        utils.check_dir(anndir)
        treedir = path.join(tree_root, str(tree_id)[:2])
        utils.check_dir(treedir)

        to_keep = []
        for n in tree.get_leaves():
            match = ens_RE.match(n.name)
            if match is not None and match.group()[:-1] in species:
                to_keep.append(n.name)

        if not len(to_keep):
            continue

        tree.prune(to_keep)

        tree.write(outfile=path.join(treedir, "{}.nh".format(tree_id)), format=5)
        tree.write(outfile=path.join(anndir, "{}.nh".format(tree_id)), features=["D", "species"], format=1)
        print "-"*40
        trees_fixed = process_tree(tree, fix=True)
        # print tree.get_ascii(show_internal=True)
        seqsets, subtrees = [], []
        for tree_fixed in trees_fixed:
            # t_seqsets, t_subtrees = split_tree_gregj(tree_fixed, tcc)
            t_seqsets, t_subtrees = split_tree(tree_fixed, tcc)
            # t_seqsets, t_subtrees = baseline_tcc(tree_fixed, tcc)
            ## Tree formatting stu
            # ts = ete3.TreeStyle()
            # ts.show_leaf_name = False
            # layout = make_layout(t_seqsets)
            # tree_fixed.show(layout=layout, tree_style=ts)
            
            for seqset in t_seqsets:
                if len(seqset) < 21:
                    sys.exit(-1)
            seqsets.extend(t_seqsets)
            subtrees.extend(t_subtrees)

        print subtrees
        outdir = path.join(out_root, args.clade, str(tree_id)[:2])
        outdir_ann = path.join(out_root_ann, args.clade, str(tree_id)[:2])
        utils.check_dir(outdir)
        utils.check_dir(outdir_ann)

        # Treevis
        # imgdir = path.join(img_root, args.clade)
        # utils.check_dir(imgdir)
        # imgfile = path.join(imgdir, "{}.pdf".format(tree_id))
        # tree.render(imgfile, layout=layout)

        set_id = 1
        for seqset, subtree in zip(seqsets, subtrees):
            # print subtree
            # taxa = [ n.name for n in subtree.get_leaves() if (ens_RE.match(n.name).group()[:-1] in species) ]

            # subtree.prune(taxa)
            # print "Pruned"
            # print subtree
            outfile = open(path.join(outdir, "{}_{}.tab".format(tree_id, set_id)), 'w')
            for seqid in seqset:
                print >>outfile, '\t'.join(seqid)

            # try:
            subtree.write(outfile=path.join(outdir, "{}_{}.nh".format(tree_id, set_id)), 
                          format=5)
            subtree.write(outfile=path.join(outdir_ann, "{}_{}.nh".format(tree_id, set_id)), 
                          format=1, features=["D", "species"])
            #except AttributeError, e:
            #    print >>sys.stderr, e

            set_id += 1

        tree_id += 1

    global abandoned
    print abandoned

if __name__ == "__main__":
    main()
