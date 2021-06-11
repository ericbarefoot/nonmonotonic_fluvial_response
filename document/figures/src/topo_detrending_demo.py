# a script to measure relief on floodplain
# Eric Barefoot
# Nov 2020

import os
import numpy as np
import pandas as pd
import h5py as h
from tqdm import tqdm

dir_path = os.getcwd()

import sys
sys.path.insert(0, dir_path)

import codebase.topo_analysis_functions as ttools
import codebase.cube.slicedice as cut

import matplotlib
from matplotlib import pyplot as plt

# matplotlib.use('Qt5Agg')

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')

topo15 = np.copy(agu_data['topoqv15'])
topo30 = np.copy(agu_data['topoqv30'])

topo15[:, 5:11, 147:160] = np.nan
topo30[:, 5:11, 147:160] = np.nan
timepick = 50

t15id = agu_data['topoqv15'].attrs['IDs'][timepick]
t30id = agu_data['topoqv30'].attrs['IDs'][timepick]
ozt15 = meta_data_19[meta_data_19.linkID == t15id].oceanZ.to_numpy()
ozt30 = meta_data_19[meta_data_19.linkID == t30id].oceanZ.to_numpy()

apex15 = (78,78)
c115 = (150,5)
c215 = (6,140)
dtheta15 = 0.05
ntheta15 = 20

tdest15 = ttools.detrend_cone(topo15[timepick], apex15, c115, c215, dtheta15)
tdest30 = ttools.detrend_cone(topo30[timepick], apex15, c115, c215, dtheta15)

f = plt.figure(figsize = (6,8))

gg = f.add_gridspec(
    ncols = 2, nrows = 4, height_ratios = [5, 1, 1, 5],
    wspace=0.15, hspace=0.5
) #  width_ratios = [1,1,1]

dd = [topo15[timepick]- ozt15/1000, tdest15]
tt = ['Scaled Topography\n$Q_v = 1.5$', 'Detrended Topography\n$Q_v = 1.5$']

base_fontsize = 8

top = f.add_subplot(gg[0, 0])
dtl = f.add_subplot(gg[0, 1])
topH = f.add_subplot(gg[1, 0])
dtlH = f.add_subplot(gg[1, 1], sharex = topH, sharey = topH)

ax = (top, dtl)

for i, dif, tit in zip(ax, dd, tt):
    # i.set_facecolor('#000000')
    i.set_title(tit, fontsize = base_fontsize, pad = -10, y = 1)
    i.spines['top'].set_visible(False)
    i.spines['bottom'].set_visible(False)
    i.spines['right'].set_visible(False)
    i.spines['left'].set_visible(False)
    i.set_xticks(np.arange(0, 2.6, 0.5))
    i.set_xticklabels(i.get_xticks(), fontsize = base_fontsize - 2)
    i.set_xlabel('Distance (m)', fontsize = base_fontsize - 1)
    i.set_yticks([])
    topomap = i.imshow(dif,
        extent = [0, (dif.shape[1] / 200), 0, (dif.shape[0] / 200)],
        origin = 'lower', cmap = 'Greys_r', vmin = -0.02, vmax = 0.04)
ax[0].set_yticks(np.arange(0, 2.6, 0.5))
ax[0].set_yticklabels(ax[0].get_yticks(), fontsize = base_fontsize - 2)
ax[0].set_ylabel('Distance (m)', fontsize = base_fontsize - 1)

k = 0
for i, j in zip((topH, dtlH), dd):
    i.hist(j.flatten(), bins = np.linspace(-0.02, 0.03, 75), color = 'k', density = True)
    i.set_xlim(-0.02, 0.03)
    i.spines['top'].set_visible(False)
    i.spines['right'].set_visible(False)
    i.set_xticks(np.arange(-0.02, 0.031, 0.01))
    i.set_xticklabels(
        np.around(i.get_xticks() * 1000, decimals = 0),
        fontsize = base_fontsize - 2
    )
    i.set_yticks(np.arange(0, 251, 50))
    i.set_yticklabels(i.get_yticks(), fontsize = base_fontsize - 2)

topH.set_xlabel('Scaled Elevation (mm)\n$Q_v = 1.5$', fontsize = base_fontsize - 1)
dtlH.set_xlabel('Detrended Elevation (mm)\n$Q_v = 1.5$', fontsize = base_fontsize - 1)
topH.set_ylabel('Frequency', fontsize = base_fontsize - 1)

dtlH.axvline(np.nanquantile(dd[1].flatten(), 0.95), ymin = 0, ymax = 0.75, color = 'r')
dtlH.axvline(np.nanquantile(dd[1].flatten(), 0.25), ymin = 0, ymax = 0.75, color = 'r')

up = np.nanquantile(dd[1].flatten(), 0.95)
lw = np.nanquantile(dd[1].flatten(), 0.25)

# print(up)

dd30 = [topo30[timepick]- ozt30/1000, tdest30]
tt30 = ['Scaled Topography\n$Q_v = 3$', 'Detrended Topography\n$Q_v = 3$']


top30 = f.add_subplot(gg[3, 0])
dtl30 = f.add_subplot(gg[3, 1])
topH30 = f.add_subplot(gg[2, 0], sharex = topH, sharey = topH)
dtlH30 = f.add_subplot(gg[2, 1], sharex = topH, sharey = topH)

ax30 = (top30, dtl30)

for i, dif, tit in zip(ax30, dd30, tt30):
    i.set_title(tit, fontsize = base_fontsize, pad = -10, y = 1)
    # i.set_facecolor('#000000')
    i.spines['top'].set_visible(False)
    i.spines['bottom'].set_visible(False)
    i.spines['right'].set_visible(False)
    i.spines['left'].set_visible(False)
    i.set_xticks(np.arange(0, 2.6, 0.5))
    i.set_xticklabels(i.get_xticks(), fontsize = base_fontsize - 2)
    i.set_xlabel('Distance (m)', fontsize = base_fontsize - 1)
    i.set_yticks([])
    topomap30 = i.imshow(dif,
        extent = [0, (dif.shape[1] / 200), 0, (dif.shape[0] / 200)],
        origin = 'lower', cmap = 'Greys_r', vmin = -0.02, vmax = 0.04)
top30.set_yticks(np.arange(0, 2.6, 0.5))
top30.set_yticklabels(ax30[0].get_yticks(), fontsize = base_fontsize - 2)
top30.set_ylabel('Distance (m)', fontsize = base_fontsize - 1)

axH30 = [topH30, dtlH30]
for i, j in zip(axH30, dd30):
    i.hist(j.flatten(), bins = np.linspace(-0.02, 0.03, 75), color = 'k', density = True)
    i.set_xlim(-0.02, 0.03)
    i.spines['top'].set_visible(False)
    i.spines['right'].set_visible(False)
    i.set_xticks(np.arange(-0.02, 0.031, 0.01))
    i.set_xticklabels(
        np.around(i.get_xticks() * 1000, decimals = 0),
        fontsize = base_fontsize - 2
    )
    # i.set_xlabel('Detrended Elevation $Q_v = 3$ (mm)', fontsize = base_fontsize - 1)
    i.set_yticks(np.arange(0, 251, 50))
    i.set_yticklabels(i.get_yticks(), fontsize = base_fontsize - 2)

topH30.set_xlabel('Scaled Elevation (mm)\n$Q_v = 3$', fontsize = base_fontsize - 1)
dtlH30.set_xlabel('Detrended Elevation (mm)\n$Q_v = 3$', fontsize = base_fontsize - 1)

dtlH30.axvline(np.nanquantile(dd30[1].flatten(), 0.95), ymin = 0, ymax = 1, color = 'r')
dtlH30.axvline(np.nanquantile(dd30[1].flatten(), 0.25), ymin = 0, ymax = 1, color = 'r')

up30 = np.nanquantile(dd30[1].flatten(), 0.95)
lw30 = np.nanquantile(dd30[1].flatten(), 0.25)

# print(up30)
# print(lw30)

topH30.set_ylabel('Frequency', fontsize = base_fontsize - 1)

cax15 = f.add_axes([0.88, 0.72, 0.01, 0.2])
cbar15 = f.colorbar(
    topomap, cax = cax15,
    ticks = np.round(np.arange(-0.02, 0.041, 0.01), 2)
)
cbar15.ax.set_yticklabels(cbar15.ax.get_yticks()*1000, fontsize = base_fontsize-2)
cbar15.ax.set_ylabel('Scaled Elevation (mm)', fontsize = base_fontsize-1)

cax30 = f.add_axes([0.88, 0.10, 0.01, 0.2])
cbar30 = f.colorbar(
    topomap30, cax = cax30,
    ticks = np.round(np.arange(-0.02, 0.041, 0.01), 2)
)
cbar30.ax.set_yticklabels(cbar30.ax.get_yticks()*1000, fontsize = base_fontsize-2)
cbar30.ax.set_ylabel('Scaled Elevation (mm)', fontsize = base_fontsize-1)

top.annotate('conical shape\ndominates roughness', xy=(0.65, 0.45),  xycoords='data',
            xytext=(1.4, 0.45), textcoords='data',
            fontsize = 8, color = 'r',
            arrowprops=dict(
                color='r', alpha = 0.75,
                arrowstyle="->",
                # connectionstyle="arc3, rad=0.5",
                relpos = (0, 0.5)
            ),
            horizontalalignment='left', verticalalignment='top',
            )

dtl.annotate('levees\nremain', xy=(0.85, 1.2),  xycoords='data',
            xytext=(0.1, 1.5), textcoords='data',
            fontsize = 8, color = 'r',
            arrowprops=dict(
                color='r', alpha = 0.75,
                arrowstyle="->",
                # connectionstyle="arc3, rad=0.5",
                relpos = (1, 0.5)
            ),
            horizontalalignment='right', verticalalignment='top',
            )


dtl30.annotate('smooth\nsurface', xy=(0.4, 1.8),  xycoords='data',
            xytext=(-0.1, 0.75), textcoords='axes fraction',
            fontsize = 8, color = 'r',
            arrowprops=dict(
                color='r', alpha = 0.75,
                arrowstyle="->",
                # connectionstyle="arc3, rad=0.5",
                relpos = (1, 0.5)
            ),
            horizontalalignment='right', verticalalignment='top',
            )


dtlH30.annotate('small\nroughness', xy=(np.mean([up30, lw30]), 250), xycoords='data',
            xytext=(0.015, 200), textcoords='data',
            fontsize = 8, color = 'r',
            arrowprops=dict(
                color='r', alpha = 0.75,
                arrowstyle="-[, widthB = 1.3, lengthB = 0.3",
                connectionstyle="arc,angleA=0,angleB=90,armA=10,armB=10,rad=5",
                relpos = (0, 0.5)
            ),
            horizontalalignment='left', verticalalignment='bottom',
            )

# print(np.mean([up, lw]))
#
dtlH.annotate('large\nroughness', xy=(np.mean([up, lw]), 200), xycoords='data',
            xytext=(0.015, 200), textcoords='data',
            fontsize = 8, color = 'r',
            arrowprops=dict(
                color='r', alpha = 0.75,
                arrowstyle="-[, widthB = 2, lengthB = 0.3",
                # ,
                connectionstyle="arc,angleA=0,angleB=90,armA=10,armB=10,rad=5",
                relpos = (0, 0.5)
            ),
            horizontalalignment='left', verticalalignment='bottom',
            )

plt.annotate('', xy=(0.05, 0.62),  xycoords='figure fraction',
            xytext=(0.05, 0.72), textcoords='figure fraction',
            fontsize = 8, color = 'r',
            arrowprops=dict(
                color='k', alpha = 0.75,
                arrowstyle="->",
                connectionstyle="arc3, rad=0.5",
                relpos = (1, 0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

plt.annotate('', xy=(0.05, 0.4),  xycoords='figure fraction',
            xytext=(0.05, 0.3), textcoords='figure fraction',
            fontsize = 8, color = 'r',
            arrowprops=dict(
                color='k', alpha = 0.75,
                arrowstyle="->",
                connectionstyle="arc3, rad=-0.5",
                relpos = (1, 0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

f.subplots_adjust(left = 0.075, bottom = 0.075, top = 0.95)

f.savefig(fname = os.path.join(dir_path, 'analysis/chapter_2/figures/output/detrending_procedure.pdf'), transparent = False, dpi = 200)

agu_data.close()
