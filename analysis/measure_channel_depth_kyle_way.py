# a script to measure channel depth and width
# Eric Barefoot
# Nov 2020

# This needs very clean channel masks to work well. At the moment, the masks
# I have are not clean enough, and end up generating totally spurious results.
# may be worth it to go through and hand correct masks for this.

import os
import numpy as np
import pandas as pd
import h5py as h
from tqdm import tqdm

import matplotlib
from matplotlib import pyplot as plt

import codebase.cube.slicedice as cut

dir_path = os.getcwd()

meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))
meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')

# import the topographic surfaces

topo10 = np.copy(agu_data['topoqv10'])
topo15 = np.copy(agu_data['topoqv15'])
topo30 = np.copy(agu_data['topoqv30'])

msk10 = np.copy(agu_data['manualChannelMasksqv10'])
msk15 = np.copy(agu_data['manualChannelMasksqv15'])
msk30 = np.copy(agu_data['manualChannelMasksqv30'])
# used to do it with the auto-generated masks. Was a mess
# msk10 = np.copy(agu_data['channelMasksqv10'])
# msk15 = np.copy(agu_data['channelMasksqv15'])
# msk30 = np.copy(agu_data['channelMasksqv30'])

topstr = ['topoqv10', 'topoqv15', 'topoqv30']
mskstr = ['manualChannelMasksqv10', 'manualChannelMasksqv15', 'manualChannelMasksqv30']
# mskstr = ['channelMasksqv10', 'channelMasksqv15', 'channelMasksqv30']
topdat = [topo10, topo15, topo30]
mskdat = [msk10, msk15, msk30]

mtches = []
j = 0
for istr, jstr, idat, jdat in zip(topstr, mskstr, topdat, mskdat):
    mtches.append(np.isin(agu_data[jstr].attrs['IDs'], agu_data[istr].attrs['IDs'] - 2))
    j += 1

apex12 = (79, 87)
apex19 = (78, 78)

qvs = [1, 1.5, 3]
metas = [meta_data_12, meta_data_19, meta_data_19]
dataIDs = [agu_data['topoqv10'].attrs['IDs'], agu_data['topoqv15'].attrs['IDs'], agu_data['topoqv30'].attrs['IDs']]
apices = [apex12, apex19, apex19]
chanDimData = {k: [] for k in ['qv', 'time', 'radius', 'width', 'depth', 'xrb', 'yrb', 'xlb', 'ylb', 'tID']}

rsweep = np.arange(0.1, 1, 0.2)

for topo, meta, apex, qv, match, mask, mstr in tqdm(zip(topdat, metas, apices, qvs, mtches, mskdat, mskstr)):
    for r in (200 * rsweep).astype('int'):
        sec, x, y = cut.circular_slice(topo, apex, r)
        nn = np.any(~np.isnan(sec), axis = 0)
        dx = np.diff(x[nn], prepend = x[nn][0])
        dy = np.diff(y[nn], prepend = y[nn][0])
        d = np.sqrt(dx * dx + dy * dy).cumsum()
        maskI = mask  # [match, :]
        mskSlc = maskI[:, x[nn], y[nn]]
        mskSlc[:, -1] = False
        mskDiff = np.diff(mskSlc.astype('int'), axis = -1, prepend = 0)
        rbs = np.where(mskDiff == 1)[1]
        rbsT = np.where(mskDiff == 1)[0]
        lbs = np.where(mskDiff == -1)[1] - 1
        lbsT = np.where(mskDiff == -1)[0]
        for rb, lb, t in zip(rbs, lbs, rbsT):
            if rb != lb:
                xrb, yrb = x[nn][rb], y[nn][rb]
                xlb, ylb = x[nn][lb], y[nn][lb]
                chanDimData['xrb'].append(xrb)
                chanDimData['yrb'].append(yrb)
                chanDimData['xlb'].append(xlb)
                chanDimData['ylb'].append(ylb)
                chanDimData['width'].append(np.sqrt((xrb - xlb)**2 + (yrb - ylb)**2) / 200)
                chanDimData['qv'].append(qv)
                chanDimData['radius'].append(r / 200)
                chanslc = np.arange(rb, lb)
                toposec = sec[t, nn]
                chanDimData['depth'].append(np.ptp(toposec[chanslc]))
                mskLink = agu_data[mstr].attrs['IDs'][t]
                # mskLink = agu_data[mstr].attrs['IDs'][match][t]
                chanDimData['tID'].append(mskLink)
                time = meta[meta.linkID == mskLink].runtime.to_numpy()
                chanDimData['time'].extend(time)

chandat = pd.DataFrame(chanDimData)

chandat.to_csv(os.path.join(dir_path, 'data/derived_data/channeldimensionsxy.csv'), index = False)

agu_data.close()
