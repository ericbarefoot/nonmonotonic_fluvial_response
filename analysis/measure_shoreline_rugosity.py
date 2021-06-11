# a script to make synstrat datasets
# Eric Barefoot
# Nov 2020

import os
import numpy as np
import pandas as pd
import h5py as h
from tqdm import tqdm

from skimage.morphology import binary_erosion as be
from matplotlib.path import Path

import codebase.cube.slicedice as cut

def make_roiMask(arr, pts):
    yy, xx = arr.shape
    x, y = np.meshgrid(np.arange(xx), np.arange(yy)) # make a canvas with coordinates
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x,y)).T
    p = Path(pts) # make a polygon
    grid = p.contains_points(points)
    mask = grid.reshape(yy,xx) # now you have a mask with points inside a polygon
    return mask

def get_shore_verts(boolArr, backwallROI, xwall, ywall):
    sh1 = be(boolArr) ^ boolArr
    sh1[xwall, :] = False
    sh1[:, ywall] = False
    bwmsk = make_roiMask(sh1, backwallROI)
    sh1[bwmsk] = False
    return np.where(sh1)

def get_rugosity(x, y, apex):
    xa, ya = apex
    xd = x - xa
    yd = y - ya
    rs = np.sqrt(xd * xd + yd * yd)
    mr = np.mean(rs)
    n = len(rs)
    return np.sqrt(np.sum(((rs - mr) / mr)**2) / n)

dir_path = os.getcwd()

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))
meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'a')

topo10 = agu_data['topoqv10']
topo15 = agu_data['topoqv15']
topo30 = agu_data['topoqv30']

backWallROI10 = [
    (159.65764557773477, 2.512843820692325),
    (160.9917367528929, 8.38284499138814),
    (8.63852454983325, 162.07014836960593),
    (2.234886909074177, 162.33696660463755),
    (2.501705144105806, 155.1328742587836),
    (151.91991676181755, 3.313298525787218)
]
apex10 = (79,87)

imax = agu_data['topoqv10'].shape[0]

tt = np.zeros(imax)
R = np.zeros(imax)
qv = [1] * imax

for i in tqdm(range(imax)):
    y, x = get_shore_verts(agu_data['topsetMasksqv10'][i*2], backWallROI10, 6, 5)
    R[i] = get_rugosity(x, y, apex10)
    ID = agu_data['topsetMasksqv10'].attrs['IDs'][i*2]
    tt[i] = meta_data_12[np.isin(meta_data_12.linkID, ID)].runtime.to_numpy()

rugData10 = {'rugosity': R, 'time': tt, 'qv': qv}
rugdat10 = pd.DataFrame(rugData10)

backWallROI15 = [
    (147.42873693244786, 2.1601899486271634),
    (159.97708031372878, 2.4390420237667456),
    (161.37134068942663, 8.852639751976994),
    (8.281551437799365, 162.77898522902296),
    (2.1468057847286843, 157.48079580137102),
    (4.098770310705717, 146.32671279578796)
]
apex15 = (78,78)

imax = agu_data['topoqv15'].shape[0]

tt = np.zeros(imax)
R = np.zeros(imax)
qv = [1.5] * imax

for i in tqdm(range(imax)):
    y, x = get_shore_verts(agu_data['topsetMasksqv15'][i*2], backWallROI15, 6, 6)
    R[i] = get_rugosity(x, y, apex15)
    ID = agu_data['topsetMasksqv15'].attrs['IDs'][i*2]
    tt[i] = meta_data_19[np.isin(meta_data_19.linkID, ID)].runtime.to_numpy()

rugData15 = {'rugosity': R, 'time': tt, 'qv': qv}
rugdat15 = pd.DataFrame(rugData15)

imax = agu_data['topoqv30'].shape[0]

tt = np.zeros(imax)
R = np.zeros(imax)
qv = [3] * imax

for i in tqdm(range(imax)):
    y, x = get_shore_verts(agu_data['topsetMasksqv30'][i*2], backWallROI15, 6, 6)
    R[i] = get_rugosity(x, y, apex15)
    ID = agu_data['topsetMasksqv30'].attrs['IDs'][i*2]
    tt[i] = meta_data_19[np.isin(meta_data_19.linkID, ID)].runtime.to_numpy()

rugData30 = {'rugosity': R, 'time': tt, 'qv': qv}
rugdat30 = pd.DataFrame(rugData30)


allRugs = rugdat10
rugdatasets_to_add = [rugdat15, rugdat30]

for i in rugdatasets_to_add:
    allRugs = allRugs.append(i)

allRugs.to_csv(os.path.join(dir_path, 'data/derived_data/shorelineRugosity.csv'), index = False)

agu_data.close()
