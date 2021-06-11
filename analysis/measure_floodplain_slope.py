# a script to measure floodplain slope for intervals in tulane data
# Eric Barefoot
# Nov 2020

import os
import numpy as np
import pandas as pd
import h5py as h
from tqdm import tqdm

import codebase.topo_analysis_functions as ttools

import matplotlib
from matplotlib import pyplot as plt

dir_path = os.getcwd()

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))
meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')

topo10 = np.copy(agu_data['topoqv10'])
topo15 = np.copy(agu_data['topoqv15'])
topo30 = np.copy(agu_data['topoqv30'])

topo15[:, 5:11, 147:160] = np.nan
topo30[:, 5:11, 147:160] = np.nan

ntheta = 200

apex10 = (79,87)
c110 = (157,6)
c210 = (5,160)

nsurf10 = topo10.shape[0]
slps10 = np.zeros((nsurf10, ntheta))
angs10 = np.zeros((nsurf10, ntheta))

for i in tqdm(range(topo10.shape[0])):
    slps10[i], angs10[i] = ttools.measure_slopes(topo10[i], apex10, c110, c210, ntheta)

apex15 = (78,78)
c115 = (150,5)
c215 = (6,140)

nsurf15 = topo15.shape[0]
slps15 = np.zeros((nsurf15, ntheta))
angs15 = np.zeros((nsurf15, ntheta))

for i in tqdm(range(topo15.shape[0])):
    slps15[i], angs15[i] = ttools.measure_slopes(topo15[i], apex15, c115, c215, ntheta)

nsurf30 = topo30.shape[0]
slps30 = np.zeros((nsurf30, ntheta))
angs30 = np.zeros((nsurf30, ntheta))

for i in tqdm(range(topo30.shape[0])):
    slps30[i], angs30[i] = ttools.measure_slopes(topo30[i], apex15, c115, c215, ntheta)

slopedata = [slps10, slps15, slps30]
angledata = [angs10, angs15, angs30]
qvs = [1, 1.5, 3]
dataIDs = [agu_data['topoqv10'].attrs['IDs'], agu_data['topoqv15'].attrs['IDs'], agu_data['topoqv30'].attrs['IDs']]
metas = [meta_data_12, meta_data_19, meta_data_19]

SlpData = {k: [] for k in ['angle', 'slope', 'qv', 'time']}

for aa, ss, theqv, ids, met in zip(angledata, slopedata, qvs, dataIDs, metas):
    ts = met[np.isin(met.linkID, ids)].runtime.to_numpy()
    for tt in range(len(ts)):
        SlpData['qv'].extend([theqv] * len(aa[tt]))
        SlpData['angle'].extend(aa[tt])
        SlpData['slope'].extend(ss[tt])
        SlpData['time'].extend([ts[tt]] * len(aa[tt]))

for i in ['angle', 'slope', 'qv', 'time']:
    print(len(SlpData[i]))

slpdat = pd.DataFrame(SlpData)

slpdat.to_csv(os.path.join(dir_path, 'data/derived_data/FPSlopes.csv'), index = False)

agu_data.close()
