# a script to measure relief on floodplain
# Eric Barefoot
# Nov 2020

import os
import numpy as np
import pandas as pd
import h5py as h
from tqdm import tqdm

import matplotlib
from matplotlib import pyplot as plt

dir_path = os.getcwd()

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))
meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'a')

# import the detrended surfaces

topo10DT = np.copy(agu_data['topoDTqv10'])
topo15DT = np.copy(agu_data['topoDTqv15'])
topo30DT = np.copy(agu_data['topoDTqv30'])
topo10 = np.copy(agu_data['topoqv10'])
topo15 = np.copy(agu_data['topoqv15'])
topo30 = np.copy(agu_data['topoqv30'])

tlen = topo10.shape[0]

rdata10 = np.zeros(tlen)
rdata15 = np.zeros(tlen)
rdata30 = np.zeros(tlen)

dd = [topo10, topo15, topo30]
ddDT = [topo10DT, topo15DT, topo30DT]
rr = [rdata10, rdata15, rdata30]
rrDT = [np.zeros(tlen), np.zeros(tlen), np.zeros(tlen)]
qvs = [1, 1.5, 3]
metas = [meta_data_12, meta_data_19, meta_data_19]
dataIDs = [agu_data['topoqv10'].attrs['IDs'], agu_data['topoqv15'].attrs['IDs'], agu_data['topoqv30'].attrs['IDs']]

Reldata = {k: [] for k in ['short_wave_relief', 'long_wave_relief', 'qv', 'time']}

for ii, idt, jj, rdt, kk, met, ids in zip(dd, ddDT, rr, rrDT, qvs, metas, dataIDs):
    for i in np.arange(tlen):
        # m = np.nanmean(ii[i])
        # ma = np.nanmax(ii[i])
        # jj[i] = ma - m
        # jj[i] = np.nanstd(ii[i])
        # jj[i] = np.ptp([~np.isnan(ii[i])])
        q1 = np.quantile(ii[i][~np.isnan(ii[i])], 0.25)
        q2 = np.quantile(ii[i][~np.isnan(ii[i])], 0.95)
        q1DT = np.quantile(idt[i][~np.isnan(idt[i])], 0.25)
        q2DT = np.quantile(idt[i][~np.isnan(idt[i])], 0.95)
        jj[i] = q2 - q1
        rdt[i] = q2DT - q1DT
    Reldata['qv'].extend([kk] * tlen)
    Reldata['long_wave_relief'].extend(jj)
    Reldata['short_wave_relief'].extend(rdt)
    Reldata['time'].extend(met[np.isin(met.linkID, ids)].runtime.to_numpy())


reldat = pd.DataFrame(Reldata)

reldat.to_csv(os.path.join(dir_path, 'data/derived_data/FPrelief.csv'), index = False)


# f, a = plt.subplots(nrows = 3,tight_layout = True, sharex = True)

# for i in a:
#     i.cla()
#
# # for i, j in zip(a,dd):
# #     oo = i.hist(j[99].flatten(), bins = 30)
#
# for i, j in zip(a,rr):
#     # oo = i.hist(j, bins = 30, range = (0.02, 0.05))
#     oo = i.hist(j, bins = 10)



# f.show()

agu_data.close()
