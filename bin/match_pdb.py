import sys
from argparse import ArgumentParser
import urllib, urllib2, httplib
import glob
from os import path

url = 'http://www.uniprot.org/uploadlists/'

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
'format':'txt',
'file': '\n'.join(all_ids)
}

print all_ids
data = urllib.urlencode(params)
request = urllib2.Request(url, data)
contact = "gregs@ebi.ac.uk" # Please set your email address here to help us debug in case of problems.
request.add_header('User-Agent', 'libwww-perl {}'.format(contact))
request.add_header('Content-Type', 'form-data')
response = urllib2.urlopen(request)

# First we fetch all Uniprot translations for the Ensembl IDs we need.
# Then go through *all* of SIFTS to find matching records.
# Finally, we output the matching ENSEMBL-PDB pairs

# Uniprot -> Ensembl
response.readline() # Skip header
for l in response:
    print l
    f = l.rstrip().split()

    if id_map.get(f[1]):
        print >>sys.stderr, f[1], "already in the map!"

    id_map[f[1]] = f[0]


# Download from the SIFTS website
sifts_file = open(args.siftsfile)
# From submission to the Uniprot translation service
out_file = open(args.outfile, 'w')

sp2pdb = {}
sifts_file.readline() # Skip header
for l in sifts_file:
    f = l.rstrip().split('\t')
    spid = f[2]

    if spid not in id_map:
        print >>sys.stderr, "Skipping", spid, l
        continue

    # New way of generating the map without making any decisions as to what the best 
    # structure is.
    sp2pdb[ id_map[spid] ] = f[:2] + [ spid ]
    print >>out_file, '\t'.join([ id_map[ spid ] ] + f[:2] + [ spid ])

# TODO Could work out which Ensembl IDs didn't get a mapping from Uniprot
# -- By set difference for example.
print >>sys.stderr, "All IDs", len(all_ids)
print >>sys.stderr, "ID map", len(id_map)
print >>sys.stderr, "SP2PDB", len(sp2pdb.keys())
print >>sys.stderr, "all IDs - SP2PDB", len(set(all_ids).difference(sp2pdb.keys()))
