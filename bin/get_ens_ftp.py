import sys
import pickle
from ftplib import FTP
from os import path
from argparse import ArgumentParser


argparser = ArgumentParser()

argparser.add_argument('--clade', metavar='clade_name', type=str, required=True)
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
argparser.add_argument('--release', metavar='version', type=str, required=True)
argparser.add_argument('--pep', action="store_true")
argparser.add_argument('--cds', action="store_true")
argparser.add_argument('--species_cache', metavar='pickle_file', type=str, required=True)

args = argparser.parse_args()

clades_pickle = args.species_cache
clades = pickle.load(open(clades_pickle))
species_list = [ sp.lower().replace(' ', '_') for sp in clades[args.clade] ]

ftp = FTP('ftp.ensembl.org')
ftp.login()
ftp.cwd('/pub/release-{}/fasta/'.format(args.release))
# species_list = ftp.nlst()

if args.pep == True:
  print >>sys.stderr, "Running in pep mode..."
  dir = '/pep'
else:
  print >>sys.stderr, "Running in cds mode..."
  dir = '/cds'

for species in species_list:
  for f in ftp.nlst(species + dir):
    if f.endswith('all.fa.gz'):
      print f
      outfile = path.join(args.outroot, path.basename(f))
      ftp.retrbinary('RETR ' + f, open(outfile, "w").write)
      break
