import os
import numpy as np
import numpy.ma as ma
import pandas as pd
import h5py as h
from tqdm import tqdm


import codebase.topo_analysis_functions as ttools
import codebase.cube.slicedice as cut

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

dir_path = os.getcwd()

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))
meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r+')

apex19 = (78, 78)
apex12 = (79, 87)

def get_completeness_slice(holes, apex, r, rLabel, qv):
    slc, xx, yy = cut.circular_slice(holes, apex, int(200 * r))
    slc = np.tile(slc, (10, 1))
    slcmsk = np.all(slc, axis = 0)
    total = (~slcmsk).sum()
    slc = ma.masked_array(slc, mask = np.tile(slcmsk, (slc.shape[0], 1)))
    compData = {k: [] for k in ['f', 'dt', 'qv', 'r', 'rLabel', 'sid']}
    for j in tqdm(range(1, slc.shape[0])):
        # j = 25  # 1 through 100
        f_tmp = np.zeros((j, slc.shape[1]))
        f_tmp = ma.masked_array(f_tmp, mask = np.tile(slcmsk, (j, 1)))
        possible_sections = np.zeros((j))
        jind, lag = (0, 0)
        # for jind, lag in enumerate(range(-(j - 1), 1, 1)):
            # lag = -20  # -(j - 1) through 0
        splitInds = [i for i in range(j + lag, slc.shape[0], j)]
        sp = np.array_split(slc, splitInds)
        lsp = [i.shape[0] for i in sp]
        tfsp = [i == j for i in lsp]
        n_splits = sum(tfsp)
        possible_sections[jind] = n_splits
        recombine = np.zeros((n_splits,) + slc.shape[1:])
        recombine = ma.masked_array(recombine, mask = np.tile(slcmsk, (n_splits, 1)))
        k = 0
        kk = 0
        for tf in tfsp:
            if tf:
                recombine[kk] = np.any(~sp[k], axis = 0)
                kk += 1
            k += 1
        f_tmp[jind, :] = recombine.sum(axis = 0)
#         if j == 45:
#             print(f_tmp.sum(axis = 0).data)
#             print(possible_sections.sum())
        f_tmp = f_tmp.sum(axis = 0) / possible_sections.sum()
        f = f_tmp.data[~f_tmp.mask]
#         if j == 45:
#             print(f)
        dt = np.ones(f.shape) * j
        rs = np.ones(f.shape) * round(r, 2)
        rint = np.ones(f.shape) * rLabel
        qvs = np.ones(f.shape) * qv
#         print(dt)
#         lags = np.ones(f.shape) * lag
        compData['dt'].extend(dt)
        compData['f'].extend(f)
        compData['qv'].extend(qvs)
        compData['r'].extend(rs)
        compData['rLabel'].extend(rint)
        compData['sid'].extend(np.arange(f.shape[0]))
#         compData['lag'].extend(lags)
    return pd.DataFrame(compData).sample(frac = 0.2)


apex19 = (78, 78)
apex12 = (79, 87)

rs = np.array(
    [np.linspace(0.3, 1, 6),
    np.linspace(0.3, 1.5, 6),
    np.linspace(0.3, 1.2, 6)]
)
# rs = np.array(
#     [np.array([0.75]),
#     np.array([0.75]),
#     np.array([0.75])]
# )

apices = [apex12, apex19, apex19]
qv = [1, 1.5, 3]
strats = [
    agu_data['synstratNanqv10'],
    agu_data['synstratNanqv15'],
    agu_data['synstratNanqv30']
]
compDat = []

for s, q, apex, rlst in zip(strats, qv, apices, rs):
    holes = np.isnan(np.copy(s))
    rlab = 1
    for r in rlst:
        compDat.append(get_completeness_slice(holes, apex, r, rlab, q))
        rlab += 1

allCompleteness = pd.concat(compDat)

allCompleteness.to_csv(os.path.join(dir_path, 'data/derived_data/strat_completeness.csv'), index = False)

agu_data.close()














###############################################################

# def get_completeness_slice(holes, apex, r, qv):
#     slc, xx, yy = cut.circular_slice(holes, apex, int(200 * r))
#     slcmsk = np.all(slc, axis = 0)
#     total = (~slcmsk).sum()
#     slc = ma.masked_array(slc, mask = np.tile(slcmsk, (slc.shape[0], 1)))
#     compData = {k: [] for k in ['f', 'dt', 'qv', 'r', 'lag']}
#     for j in range(1, slc.shape[0]):
#         # j = 25  # 1 through 100
#         for lag in range(-(j - 1), 1):
#             # lag = -20  # -(j - 1) through 0
#             splitInds = [i for i in range(j + lag, slc.shape[0], j)]
#             sp = np.array_split(slc, splitInds)
#             lsp = [i.shape[0] for i in sp]
#             tfsp = [i == j for i in lsp]
#             n_splits = sum(tfsp)
#             recombine = np.zeros((n_splits,) + slc.shape[1:])
#             recombine = ma.masked_array(recombine, mask = np.tile(slcmsk, (n_splits, 1)))
#             k = 0
#             kk = 0
#             for tf in tfsp:
#                 if tf:
#                     recombine[kk] = np.any(~sp[k], axis = 0)
#                     kk += 1
#                 k += 1
#             f_tmp = recombine.sum(axis = 0) / n_splits
#             f = f_tmp.data[~f_tmp.mask]
#             dt = np.ones(f.shape) * j
#             rs = np.ones(f.shape) * r
#             qvs = np.ones(f.shape) * qv
#             lags = np.ones(f.shape) * lag
#             compData['dt'].extend(dt)
#             compData['f'].extend(f)
#             compData['qv'].extend(qvs)
#             compData['r'].extend(rs)
#             compData['lag'].extend(lags)
#     return pd.DataFrame(compData)
