# a script to measure length of bankline and area of channel.
# Eric Barefoot
# Nov 2020

import os
import numpy as np
import pandas as pd
import h5py as h
from tqdm import tqdm

from skimage.util import img_as_bool
from skimage.morphology import binary_dilation as bd
from matplotlib.path import Path

# import codebase.cube.slicedice as cut

# import matplotlib
# from matplotlib import pyplot as plt
# from matplotlib.collections import LineCollection
# matplotlib.use('Qt5Agg')
# plt.ion()


def get_bank_img(boolArr, backwallROI, xwall, ywall):
    sh1 = boolArr ^ bd(boolArr)
    sh1[xwall, :] = False
    sh1[:, ywall] = False
    bwmsk = make_roiMask(sh1, backwallROI)
    sh1[bwmsk] = False
    return sh1


def make_roiMask(arr, pts):
    yy, xx = arr.shape
    x, y = np.meshgrid(np.arange(xx), np.arange(yy)) # make a canvas with coordinates
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x,y)).T
    p = Path(pts) # make a polygon
    grid = p.contains_points(points)
    mask = grid.reshape(yy, xx) # now you have a mask with points inside a polygon
    return mask


dir_path = os.getcwd()

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))
meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')

backWallROI10 = [
    (159.65764557773477, 2.512843820692325),
    (160.9917367528929, 8.38284499138814),
    (8.63852454983325, 162.07014836960593),
    (2.234886909074177, 162.33696660463755),
    (2.501705144105806, 155.1328742587836),
    (151.91991676181755, 3.313298525787218)
]
apex10 = (79,87)

imax = agu_data['manualChannelMasksqv10'].shape[0]

tt = np.zeros(imax)
P = np.zeros(imax)
A = np.zeros(imax)
qv = [1] * imax

for i in tqdm(range(imax)):
    bwmsk = img_as_bool(agu_data['manualChannelMasksqv10'][i])
    A[i] = bwmsk.sum() * 25
    bimg = get_bank_img(bwmsk, backWallROI10, 6, 5)
    P[i] = bimg.sum() * 5
    ID = agu_data['manualChannelMasksqv10'].attrs['IDs'][i]
    tt[i] = meta_data_12[np.isin(meta_data_12.linkID, ID)].runtime.to_numpy()

perimData10 = {'area_channel': A, 'channel_bank_length': P, 'time': tt, 'qv': qv}
perimdat10 = pd.DataFrame(perimData10)

backWallROI15 = [
    (147.42873693244786, 2.1601899486271634),
    (159.97708031372878, 2.4390420237667456),
    (161.37134068942663, 8.852639751976994),
    (8.281551437799365, 162.77898522902296),
    (2.1468057847286843, 157.48079580137102),
    (4.098770310705717, 146.32671279578796)
]
apex15 = (78,78)

imax = agu_data['manualChannelMasksqv15'].shape[0]

tt = np.zeros(imax)
P = np.zeros(imax)
A = np.zeros(imax)
qv = [1.5] * imax

for i in tqdm(range(imax)):
    bwmsk = img_as_bool(agu_data['manualChannelMasksqv15'][i])
    A[i] = bwmsk.sum() * 25
    bimg = get_bank_img(bwmsk, backWallROI15, 6, 6)
    P[i] = bimg.sum() * 5
    ID = agu_data['manualChannelMasksqv15'].attrs['IDs'][i]
    tt[i] = meta_data_19[np.isin(meta_data_19.linkID, ID)].runtime.to_numpy()

perimData15 = {'area_channel': A, 'channel_bank_length': P, 'time': tt, 'qv': qv}
perimdat15 = pd.DataFrame(perimData15)

imax = agu_data['manualChannelMasksqv30'].shape[0]

tt = np.zeros(imax)
P = np.zeros(imax)
A = np.zeros(imax)
qv = [3] * imax

for i in tqdm(range(imax)):
    bwmsk = img_as_bool(agu_data['manualChannelMasksqv30'][i])
    A[i] = bwmsk.sum() * 25
    bimg = get_bank_img(bwmsk, backWallROI15, 6, 6)
    P[i] = bimg.sum() * 5
    ID = agu_data['manualChannelMasksqv30'].attrs['IDs'][i]
    tt[i] = meta_data_19[np.isin(meta_data_19.linkID, ID)].runtime.to_numpy()

perimData30 = {'area_channel': A, 'channel_bank_length': P, 'time': tt, 'qv': qv}
perimdat30 = pd.DataFrame(perimData30)


allPerims = perimdat10
perimdatasets_to_add = [perimdat15, perimdat30]

for i in perimdatasets_to_add:
    allPerims = allPerims.append(i)

allPerims.to_csv(os.path.join(dir_path, 'data/derived_data/channel_areas_perimeters.csv'), index = False)

agu_data.close()
