from ftplib import FTP
from os import path

ftp = FTP('ftp.ensembl.org')
ftp.login()
ftp.cwd('/pub/release-73/fasta/')
species_list = ftp.nlst()

for species in species_list:
  if species == "ancestral_alleles":
    continue
  for f in ftp.nlst(species + '/pep'):
    if f.endswith('all.fa.gz'):
      print f
      ftp.retrbinary('RETR ' + f, open(path.basename(f), "w").write)
      break
