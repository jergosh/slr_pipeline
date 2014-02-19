import os
from os import path
import sys
import pickle
import operator
import httplib2
import urllib
import json
import time
from pprint import pprint

import ete2
import emf

http = httplib2.Http(".cache")

clade_dir = "data/clades"
img_dir = "data/img"
tree_root = "data/ens/73/trees"
clades_pickle = "data/clades.pk"

clade_names = [ "Eutheria" ]

server = "http://beta.rest.ensembl.org/"
def ens_get(ext, *args, **kwargs):
    if len(args):
        ext += '&'.join([ urllib.quote(a) for a in args]) 
        
    if len(kwargs):
        ext += '&' + urllib.urlencode(kwargs)

    print ext

    resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    
    if not resp.status == 200:
        exc = IOError(resp)
        exc.errno = resp
        raise exc

    decoded = json.loads(content)
    return decoded

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
            sp_info = ens_get("taxonomy/id/", parent)
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
    def __init__(self, clade, condition, thr):
        self.clade = clade
        self.condition = condition
        self.thr = thr

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

def split_tree(intree, tcclist):
    seqsets, subtrees = [], []

    trees = []
    to_remove = []
    for node in intree.traverse("postorder"):
        if node.dist > 100:
            node.detach()
            if len(node.get_descendants()) > 0:
                trees.append(node)
        # TODO: Does this do anything?
        elif node.is_leaf() and not getattr(node, 'species', None):
            print "REMOVING b/c of lack of species ann."
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

    for tree in trees:
        for node in tree.traverse("postorder"):
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
                if tcclist.check(node):
                    node.add_features(split=True)
                else:
                    node.add_features(split=False)

                node.add_features(done=False)

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
    all_species = ens_get("/info/species/")["species"]
    all_species_names = [ it["name"].replace("_", " ") for it in all_species ]
    all_species_names.remove("Ancestral sequences")

    if path.exists(clades_pickle):
        Clades = pickle.load(open(clades_pickle, 'rb'))
    else:
        Clades = filter_clades(all_species_names,
                               [ "Eutheria", "Glires", "Laurasiatheria", "Sauria", "Mammalia" ])
        pickle.dump(Clades, open(clades_pickle, 'wb'))

    pprint(Clades)

    TL = TCCList()
    TL.add(TCC(Clades["Eutheria"], operator.ge, 0.6))

    clade = "Eutheria"
    outroot = path.join("data/ens/73/seqsets/", clade)
    if not path.exists(outroot):
        os.mkdir(outroot)

    tree_id = 1
    for tree in emf.EMF("data/Compara.73.protein.nhx.emf"):
    # for tree in emf.EMF("/Users/greg/Downloads/Compara.nhx_trees.57.emf"):
        print tree_id
        treedir = path.join(tree_root, str(tree_id)[:2])
        if not path.exists(treedir):
            os.mkdir(treedir)

        tree.write(outfile=path.join(treedir, "{}.nh".format(tree_id)))

        seqsets, subtrees = split_tree(tree, TL)
        outdir = path.join(outroot, str(tree_id)[:2])
        if not path.exists(outdir):
            os.mkdir(outdir)

        # Treevis
        # layout = make_layout(seqsets)
        # imgfile = path.join(img_dir, "{}.pdf".format(tree_id))
        # tree.render(imgfile, layout=layout)

        set_id = 1
        for seqset, subtree in zip(seqsets, subtrees):
            outfile = open(path.join(outdir, "{}_{}.tab".format(tree_id, set_id)), 'w')
            for seqid in seqset:
                print >>outfile, '\t'.join(seqid)

            subtree.write(outfile=path.join(outdir, "{}_{}.nh".format(tree_id, set_id)))
            set_id += 1

        tree_id += 1

if __name__ == "__main__":
    main()
