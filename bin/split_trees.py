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
from Bio import Entrez
Entrez.email = "gregs@ebi.ac.uk"

http = httplib2.Http(".cache")

clade_dir = "data/clades"

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

def split_tree(tree, tcclist):
    seqsets = []
    for node in tree.traverse("postorder"):
        children = node.get_children()
        if any([ n.done for n in children ]):
            # if not all([ n.done for n in children ]):
            #     print "Warning!"
            node.add_features(done=True)
            continue

        if len(children) and all([ n.split for n in children ]):
            print "Match", node.name
            seqset = []
            for ch in children:
                seqsets.append([ (n.name, n.species) for n in ch.get_leaves() ])
            node.add_features(done=True)
        else:
            if tcclist.check(node):
                node.add_features(split=True)
            else:
                node.add_features(split=False)

            node.add_features(done=False)

    return seqsets

clades_pickle = "data/clades.pk"
def main():
    all_species = ens_get("/info/species/")["species"]
    all_species_names = [ it["name"].replace("_", " ") for it in all_species ]
    all_species_names.remove("Ancestral sequences")

    if path.exists(clades_pickle):
        Clades = pickle.load(open(clades_pickle, 'rb'))
    else:
        Clades = filter_clades(all_species_names,
                               [ "Eutheria", "Glires", "Laurasiatheria", "Sauria" ])
        pickle.dump(Clades, open(clades_pickle, 'wb'))

    pprint(Clades)

    TL = TCCList()
    TL.add(TCC(Clades["Eutheria"], operator.ge, 0.6))

    for tree in emf.EMF("/Users/greg/Downloads/Compara.73.protein.nhx.emf"):
        seqsets = split_tree(tree, TL)

if __name__ == "__main__":
    main()
