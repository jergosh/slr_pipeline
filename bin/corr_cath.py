import sys
import pandas
import argparse

import numpy as np
from scipy.stats import spearmanr
from statsmodels.sandbox.stats.multicomp import multipletests as p_adjust


def cors(df):
    nrow = len(set(df.stable_id))
    print df.cath_id.iloc[0], nrow
    if nrow < 2:
        return

    cath_matrix = df.omega.reshape(nrow, -1)

    all_ids = list(df.stable_id.iloc[range(0, df.shape[0], cath_matrix.shape[1])])
    all_ids_2 = []
    for i in df.stable_id:
        if i not in all_ids_2:
            all_ids_2.append(i)
    assert np.array_equal(all_ids, all_ids_2)

    assert set(all_ids) == set(df.stable_id)
    idx1, idx2 = [ list(i) for i in np.triu_indices(nrow, k=1) ]

    out_df = pandas.DataFrame(columns=( 'id_1', 'id_2', 'n_obs', 'cor', 'pval' ))
    for x, y in zip(idx1, idx2):
        complete_obs = [ False if (np.isnan(i) or np.isnan(j)) else True for i, j in zip(cath_matrix[x, ...], cath_matrix[y, ...])  ] 
        complete_obs = np.array(complete_obs, dtype=bool)
        
        # print all_ids[x], all_ids[y], np.sum(complete_obs)
        if np.sum(complete_obs) < 5:
            continue

        if np.std(cath_matrix[x, complete_obs]) == 0.0 or \
           np.std(cath_matrix[y, complete_obs]) == 0.0:
            print "Skipping", all_ids[x], all_ids[y]
            continue

        cor, pval = spearmanr(cath_matrix[x, complete_obs],
                              cath_matrix[y, complete_obs])
        out_df.loc[out_df.shape[0]] = [ all_ids[x],
                                        all_ids[y],
                                        np.sum(complete_obs),
                                        cor,
                                        pval ] 
    return out_df


argparser = argparse.ArgumentParser()
args = argparser.parse_args()

cath_all = pandas.read_table(open("data/cath_new.tab"), sep="\t")
cath_all.sort(columns=["cath_id"])
cath_cors = cath_all.groupby("cath_id").apply(cors)
print cath_cors
p_adjusted = pandas.Series(p_adjust(np.array(cath_cors.pval),
                                    method='fdr_bh')[1], index=cath_cors.index)
cath_cors.loc[:, 'p.adjusted'] = p_adjusted
cath_cors.to_csv("cath_cors_py_2.tab", sep="\t", na_rep="NA")
