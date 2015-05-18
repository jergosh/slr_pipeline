import sys
from os import path
import re
from argparse import ArgumentParser
from Bio import SeqIO, SeqUtils
from Bio import PDB
from Bio.PDB.DSSP import DSSP
from Bio.PDB.PDBExceptions import PDBException

import urllib
import urllib2
import StringIO
import gzip
import utils

from lxml import etree
from collections import defaultdict
from ConfigParser import _Chainmap as ChainMap

re_uniprotid = re.compile(".*\|(.*)\|.*")

# TODO move SIFTS-related stuff to a separate file?
# Caches files in 
def get_sifts(pdb_id, cache_dir):
    sifts_url = "http://www.rcsb.org/pdb/files/{}.sifts.xml.gz"

    sifts_fh = path.join(cache_dir, pdb_id+'.sifts.xml')
    if path.exists(sifts_fh):
        sifts = open(sifts_fh).read()
    else:
        try:
            request = urllib2.Request(sifts_url.format(pdb_id))
            request.add_header('User-Agent', 'Python')
            response = urllib2.urlopen(request)
        except urllib2.HTTPError, e:
            return None

        sifts_gz = StringIO.StringIO(response.read())
        sifts = gzip.GzipFile(fileobj=sifts_gz).read()
        
        utils.check_dir(cache_dir)
        with open(sifts_fh, 'w') as outfile:
            outfile.write(sifts)

    return sifts
    
    
# TODO ALl chains and IDs?
# OR: coords set up by a particular line in the aggregate SIFTS file?
# Former is probably faster
# Could return a nested dictionary with Uniprot -- PDB -- chain hierarchy
# -> ordereddict for the residue mapping
def parse_sifts(sifts, target_id):
    ns = {"ns": "http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd"}
    tree = etree.fromstring(sifts)
    ents = tree.xpath("//ns:entry/ns:entity",
                      namespaces=ns)

    result = []
    for ent in ents:
        segments = ent.xpath("//ns:segment", namespaces=ns)
        for segment in segments:
            segment_map = {}

            pdb_info = segment.xpath(".//ns:listMapRegion/ns:mapRegion/ns:db[@dbSource = 'PDB']",
                                     namespaces=ns)
            assert len(pdb_info) == 1
            pdb_id = pdb_info[0].attrib['dbAccessionId']
            pdb_chain = pdb_info[0].attrib['dbChainId']
            print "Processing", pdb_id, pdb_chain,

            uniprot_info = segment.xpath(".//ns:listMapRegion/ns:mapRegion/ns:db[@dbSource = 'UniProt']",
                                         namespaces=ns)
            if len(uniprot_info) != 1:
                print "Problem with UniProt info, skipping (" + str(uniprot_info) + ")"
                continue
            else:
                print

            uniprot_id = uniprot_info[0].attrib['dbAccessionId']
            if uniprot_id != target_id:
                continue

            for res in segment.xpath(".//ns:listResidue/ns:residue", namespaces=ns):
                res_map = res.xpath(".//ns:crossRefDb[@dbSource = 'PDB' or @dbSource = 'UniProt']",
                                    namespaces=ns)
                assert len(res_map) == 2
                assert res_map[0].attrib['dbSource'] == 'PDB'
                assert res_map[1].attrib['dbSource'] == 'UniProt'
                
                segment_map[int(res_map[1].attrib['dbResNum'])] = [pdb_id, pdb_chain, res_map[0].attrib['dbResNum']]
                # print SeqUtils.seq1(res_map[0].attrib['dbResName']), res_map[1].attrib['dbResName']
                
            result.append(segment_map)

        return result
        # for seg in segments:
        #      pdb_seq3 = ''.join([ r.attrib['dbResName'] for r in seg.xpath("//ns:listResidue/ns:residue/ns:crossRefDb[@dbSource = 'PDB']", namespaces=ns) ])
        #     pdb_seq = SeqUtils.seq1(pdb_seq3)
        #     uniprot_seq = ''.join([ r.attrib['dbResName'] for r in seg.xpath("//ns:listResidue/ns:residue/ns:crossRefDb[@dbSource = 'UniProt']", namespaces=ns) ])
 
def map_pdb():
    pass

argparser = ArgumentParser()

argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--pdbmap', metavar='pdb_map', type=str, required=True)
argparser.add_argument('--siftsdir', metavar='sifts_dir', type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

if __name__ == "__main__":
    args = argparser.parse_args()

    pdb_map_file = open(args.pdbmap)
    uniprot2pdb = defaultdict(list)

    # for l in pdb_map_file:
    for l in list(pdb_map_file):
        f = l.rstrip().split('\t')
        ens, pdb_name, chain_name, uniprot = f[:4]
        print uniprot, pdb_name, chain_name
        sifts = get_sifts(pdb_name, args.siftsdir)
        if sifts is not None:
            uniprot2pdb[(ens, uniprot)].extend(parse_sifts(sifts, uniprot))

    outfile = open(args.outfile, 'w')

    for ids, segments in uniprot2pdb.items():
        ens, uniprot = ids
        pdb_map = {}
        segments.sort(key=len, reverse=True)
        for segment in segments:
            for res in segment.items():
                if res[0] not in pdb_map:
                    pdb_map[res[0]] = res[1]

        for uniprot_coord in sorted(pdb_map):
            pdb_info = pdb_map[uniprot_coord]
            pdb_id, pdb_chain, res_id = pdb_info

            print >>outfile, '\t'.join([ ens, uniprot, str(uniprot_coord) ] + pdb_info)

