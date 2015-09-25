import os
from os import path
import httplib2
import urllib
import json
from itertools import izip_longest

import re
import pprint

from Bio import PDB

def check_dir(dirname):
    if not path.exists(dirname):
        os.mkdir(dirname)

# Ensembl REST API stuff
http = httplib2.Http(".cache")

server = "http://rest.ensembl.org/"
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


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def parse_chain(pdb_chain, pdb_begin_id, pdb_begin_ins, pdb_end_id, pdb_end_ins):
    pdb_chain = list(pdb_chain)

    for i, r in enumerate(pdb_chain):
        if r.id[1] == pdb_begin_id and r.id[2] == pdb_begin_ins:
            found_begin_id = r.id[1]
            found_begin_i = i
            break
        elif r.id[1] > pdb_begin_id:
            found_begin_id = r.id[1]
            found_begin_i = i
            break
    else:
        raise IndexError("Residue not found")

    for i, r in reversed(list(enumerate(pdb_chain))):
        if r.id[1:] == (pdb_end_id, pdb_end_ins):
            found_end_id = r.id[1]
            found_end_i = i
            break
        # We need to check if the residue is an aa because some chains contain misnumbered HETATMs 
        elif r.id[1] < pdb_end_id and r.id[0] == ' ' and PDB.is_aa(r):
            found_end_id = r.id[1]
            found_end_i = i
            break
    else:
        raise IndexError("Residue not found")

    return found_begin_id, found_begin_i, found_begin_id, found_end_i

ens_re = re.compile("^([A-Z]+)")
species_file = open("species_codes.txt")
species_map = {}

for l in species_file:
    f = l.rstrip().split('\t')
    species_map[f[0]] = f

# pprint.pprint(species_map)
def ens2species(ens_id, long_name=False):
    prefix = ens_re.match(ens_id)
    if prefix == None:
        return None

    prefix = prefix[:-1]
    if prefix in species_map:
        if long_name:
            return species_map[prefix][1]
        else:
            return species_map[prefix][2]
    else:
        return None
    
def parse_stk(fh):
    seqs = {}
    for l in fh:
        l = l.rstrip()
        if l.startswith("#"):
            continue

        if not len(l):
            continue 

        if l == "//":
            break

        parts = [x.strip() for x in l.split(" ", 1)]
        id, seq = parts

        if id not in seqs:
            seqs[id] = seq
        else:
            seqs[id] += seq

    return seqs
