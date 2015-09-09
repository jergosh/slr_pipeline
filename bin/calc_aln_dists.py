import glob
import argparse
from os import path

from Bio import AlignIO
import pandas

import numpy as np

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SubsMat.MatrixInfo import pam250

AAs = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'W', 'V', 'Y' ]

subst_matrix = np.zeros((20, 20))
for x in range(20):
    for y in range(20):
        if (AAs[x], AAs[y]) in pam250:
            subst_matrix[x, y] = pam250[(AAs[x], AAs[y])]
        else:
            subst_matrix[x, y] = pam250[(AAs[y], AAs[x])]

def site_freqs(aln):
    site_freqs = []
    for i in range(aln.get_alignment_length()):
        freqs = ProteinAnalysis(''.join(aln[:, i])).get_amino_acids_percent()
        freqs = [ freqs[aa] for aa in AAs ]
        site_freqs.append(freqs)
        
    return np.array(site_freqs)

def process_family(df, domaindir, dataset_map):
    # FIXME How to make sure we don't have to reorder the pairs to identify them
    # perhaps reorder post-hoc in R
    res_cath_id, res_id1, res_id2, res_dist = [], [], [], []
    cath_id = df.cath_id.iloc[0]
    profile_dict = {}
    print cath_id

    for i, l in df.iterrows():
        if l.stable_id not in dataset_map:
            continue
        dataset = dataset_map[l.stable_id]

        aln_file = path.join(domaindir,
                             dataset.partition('_')[0][:2],
                             l.stable_id+'_'+dataset+'.fa')
        aln = AlignIO.read(aln_file, 'fasta')
        freqs = site_freqs(aln)
        profile_dict[l.stable_id] = freqs
        
    all_ids = profile_dict.keys()

    if len(all_ids) < 2:
        return

    for id_1, id_2 in zip(*[ list(i) for i in np.triu_indices(len(all_ids), 1) ]):
        aln_1 = np.transpose(profile_dict[all_ids[id_1]])
        aln_2 = profile_dict[all_ids[id_2]]
        if len(aln_1.shape) == 1 or len(aln_2.shape) == 1:
            continue

        dist = 0.0
        for i in range(aln_1.shape[1]):
            dist += np.dot(np.dot(aln_1[:, i], subst_matrix), aln_2[i, :])

        dist /= aln_1.shape[1]
        res_cath_id.append(cath_id)
        res_id1.append(all_ids[id_1])
        res_id2.append(all_ids[id_2])
        res_dist.append(dist)
        
    return pandas.DataFrame(data={ "id_1": res_id1,
                                   "id_2": res_id2,
                                   "dist": res_dist,
                                   "cath_id": res_cath_id })

argparser = argparse.ArgumentParser()

argparser.add_argument('--domaindir', metavar="dir", type=str, required=True)
argparser.add_argument('--dataset_map', metavar="file", type=str, required=True)
argparser.add_argument('--cath_map', metavar="file", type=str, required=True)
argparser.add_argument('--outfile', metavar="output_file", type=str, required=True)


def main():
    args = argparser.parse_args()

    dataset_map = {}
    for l in open(args.dataset_map):
        f = l.rstrip().split('\t')
        dataset_map[f[0]] = f[1]
  
    cath_map = pandas.read_table(args.cath_map, sep='\t',
                                 names=["stable_id", "coords", "cath_id", "pdb_id"])
    # For each protein family, iterate over all pairs
    all_dists = cath_map.groupby('cath_id').apply(process_family, args.domaindir, dataset_map)
    all_dists.to_csv(args.outfile, sep='\t', quoting=False, index=False)


if __name__ == "__main__":
    main()