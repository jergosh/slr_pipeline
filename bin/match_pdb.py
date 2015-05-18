import sys
from argparse import ArgumentParser
import urllib, urllib2
import glob
from os import path

url = 'http://www.uniprot.org/mapping/'

argparser = ArgumentParser()

argparser.add_argument('--indir', metavar='input_dir', type=str, required=True)
argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--siftsfile', metavar='sifts_file', type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

args = argparser.parse_args()

# Glob for all Ensembl IDs in a given directory
all_ids = [ path.basename(f).partition('_')[0] for
            f in glob.glob(path.join(args.indir, args.clade, '*', 'ENS*')) ]
id_map = {}

params = {
'from':'ENSEMBL_PRO_ID',
'to':'ACC',
'format':'tab',
'query': ' '.join(all_ids)
}

data = urllib.urlencode(params)
request = urllib2.Request(url, data)
contact = "gregs@ebi.ac.uk" # Please set your email address here to help us debug in case of problems.
request.add_header('User-Agent', 'Python %s' % contact)
response = urllib2.urlopen(request)
response.readline() # Skip header
for l in response:
    f = l.rstrip().split()
    id_map[f[1]] = f[0]


# Download from the SIFTS website
sifts_file = open(args.siftsfile)
# FIXME what did I mean by this?
# From submission to the Uniprot translation service
out_file = open(args.outfile, 'w')

sp2pdb = {}
sifts_file.readline() # Skip header
for l in sifts_file:
    f = l.rstrip().split('\t')
    if f[2] not in id_map:
        # print "Skipping", f[2]
        continue

    # if f[2] in sp2pdb:
    #     pdb_len = int(f[4]) - int(f[3])
    #     pdb_len_old = int(sp2pdb[f[2]][4]) - int(sp2pdb[f[2]][3])
    #     if pdb_len > pdb_len_old:
    #         print "Replacing", sp2pdb[f[2]][0], "with", 
    #         sp2pdb[ id_map[f[2]] ] = f
    # else:
    #     sp2pdb[ id_map[f[2]] ] = f

# for i in sp2pdb.items():
#     print >>out_file, '\t'.join([ i[0] ] + i[1])

    # New way of generating the map without making any decisions as to what the best 
    # structure is.
    sp2pdb[ id_map[f[2]] ] = f
    print >>out_file, '\t'.join([ id_map[ f[2]] ] + f)

# TODO Could work out which Ensembl IDs didn't get a mapping from Uniprot
# -- By set difference for example.
print "All IDs", len(all_ids)
print "ID map", len(id_map)
print "SP2PDB", len(sp2pdb.keys())
print "all IDs - SP2PDB", len(set(all_ids).difference(sp2pdb.keys()))
