import sys
import pandas
import argparse

import numpy as np
from scipy.stats import spearmanr
from statsmodels.sandbox.stats.multicomp import multipletests as p_adjust


def cors(df):
    ncol = len(set(df.stable_id))
    print df.cath_id.iloc[0], ncol
    if ncol <= 2:
        return

    cath_matrix = df.omega.reshape(-1, ncol)
    all_ids = list(df.stable_id.iloc[range(0, df.shape[0], cath_matrix.shape[0])])
    assert set(all_ids) == set(df.stable_id)
    cors, pvals = spearmanr(cath_matrix)
    idx1, idx2 = [ list(i) for i in np.triu_indices(ncol, k=1) ]

    print cors.shape, cors[(0, 0)], pvals[(0, 0)]
    out_df = pandas.DataFrame(columns=( 'id_1', 'id_2', 'cor', 'pval' ))
    for x, y in zip(idx1, idx2):
        out_df.loc[out_df.shape[0]] = [ all_ids[x],
                                        all_ids[y],
                                        cors[(x, y)],
                                        pvals[(x, y)] ] 

    p_adjusted = pandas.Series(p_adjust(np.array(out_df.pval),
                                        method='fdr_bh')[1], index=out_df.index)
    out_df.loc[:, 'p.adjusted'] = p_adjusted

    return out_df


argparser = argparse.ArgumentParser()
args = argparser.parse_args()

cath_all = pandas.read_table(open("data/cath_new.tab"), sep="\t")
cath_all.sort(columns=["cath_id"])
cath_all.groupby("cath_id").apply(cors).to_csv("cath_cors_py.tab", sep="\t")
