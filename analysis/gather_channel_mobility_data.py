# a script to gather channel mobility data from experimental datasets.
# Eric Barefoot
# Nov 2020

import os
import numpy as np
import pandas as pd
import h5py as h
from skimage import transform as tr
from tqdm import tqdm
from codebase.plan import mobility as lm

dir_path = os.getcwd()

# print(dir_path)

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))
meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))

chan10 = agu_data['manualChannelMasksqv10']
chan15 = agu_data['manualChannelMasksqv15']
chan30 = agu_data['manualChannelMasksqv30']

# calculate channel overlap for no flooding

# print(chan10.shape)

t0, t1, ov = lm.many_overlap(chan10)

T0 = meta_data_12[np.isin(meta_data_12.linkID, chan10.attrs['IDs'])].runtime.to_numpy()[t0.astype(int)]
T1 = meta_data_12[np.isin(meta_data_12.linkID, chan10.attrs['IDs'])].runtime.to_numpy()[t1.astype(int)]

ovdat10 = pd.DataFrame({'t0': T0, 't1': T1, 'ov': ov, 'qv': 1.0})

ovdat10.to_csv(os.path.join(dir_path, 'data/derived_data/qv10overlap.csv'), index = False)

# calculate channel overlap for low intensity flooding

t0, t1, ov = lm.many_overlap(chan15)

T0 = meta_data_19[np.isin(meta_data_19.linkID, chan15.attrs['IDs'])].runtime.to_numpy()[t0.astype(int)]
T1 = meta_data_19[np.isin(meta_data_19.linkID, chan15.attrs['IDs'])].runtime.to_numpy()[t1.astype(int)]

ovdat15 = pd.DataFrame({'t0': T0, 't1': T1, 'ov': ov, 'qv': 1.5})

ovdat15.to_csv(os.path.join(dir_path, 'data/derived_data/qv15overlap.csv'), index = False)

# calculate channel overlap for high intensity flooding

t0, t1, ov = lm.many_overlap(chan30)

T0 = meta_data_19[np.isin(meta_data_19.linkID, chan30.attrs['IDs'])].runtime.to_numpy()[t0.astype(int)]
T1 = meta_data_19[np.isin(meta_data_19.linkID, chan30.attrs['IDs'])].runtime.to_numpy()[t1.astype(int)]

ovdat30 = pd.DataFrame({'t0': T0, 't1': T1, 'ov': ov, 'qv': 3.0})

ovdat30.to_csv(os.path.join(dir_path, 'data/derived_data/qv30overlap.csv'), index = False)

# calculate fluvial reworking for no flooding case

RWdata = {k: [] for k in ['t0', 't1', 'rw', 'qv']}

iterl = chan10.shape[0]
IDtimes = meta_data_12[np.isin(meta_data_12.linkID, chan10.attrs['IDs'])].runtime.to_numpy()

for ii in tqdm(range(iterl), desc = 'reworking; Qv = 1'):
    rw = 1 - lm.reworking(chan10[ii:])
    tw = np.arange(ii, iterl)
    T0 = IDtimes[ii]
    TW = IDtimes[tw]
    RWdata['qv'].extend([1] * len(rw))
    RWdata['t0'].extend([T0] * len(rw))
    RWdata['t1'].extend(TW)
    RWdata['rw'].extend(rw)

rwdat10 = pd.DataFrame(RWdata)

rwdat10.to_csv(os.path.join(dir_path, 'data/derived_data/qv10reworking.csv'), index = False)

# calculate fluvial reworking for low intenstiy flooding case

RWdata = {k: [] for k in ['t0', 't1', 'rw', 'qv']}

iterl = chan15.shape[0]
IDtimes = meta_data_19[np.isin(meta_data_19.linkID, chan15.attrs['IDs'])].runtime.to_numpy()

for ii in tqdm(range(iterl), desc = 'reworking; Qv = 1.5'):
    rw = 1 - lm.reworking(chan15[ii:])
    tw = np.arange(ii, iterl)
    T0 = IDtimes[ii]
    TW = IDtimes[tw]
    RWdata['qv'].extend([1.5] * len(rw))
    RWdata['t0'].extend([T0] * len(rw))
    RWdata['t1'].extend(TW)
    RWdata['rw'].extend(rw)

rwdat15 = pd.DataFrame(RWdata)

rwdat15.to_csv(os.path.join(dir_path, 'data/derived_data/qv15reworking.csv'), index = False)

# calculate fluvial reworking for high intenstiy flooding case

RWdata = {k: [] for k in ['t0', 't1', 'rw', 'qv']}

iterl = chan30.shape[0]
IDtimes = meta_data_19[np.isin(meta_data_19.linkID, chan30.attrs['IDs'])].runtime.to_numpy()

for ii in tqdm(range(iterl), desc = 'reworking; Qv = 3'):
    rw = 1 - lm.reworking(chan30[ii:])
    tw = np.arange(ii, iterl)
    T0 = IDtimes[ii]
    TW = IDtimes[tw]
    RWdata['qv'].extend([3] * len(rw))
    RWdata['t0'].extend([T0] * len(rw))
    RWdata['t1'].extend(TW)
    RWdata['rw'].extend(rw)

rwdat30 = pd.DataFrame(RWdata)

rwdat30.to_csv(os.path.join(dir_path, 'data/derived_data/qv30reworking.csv'), index = False)

agu_data.close()
