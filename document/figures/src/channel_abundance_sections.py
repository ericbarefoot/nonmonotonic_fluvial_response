import os
import numpy as np
import numpy.ma as ma
import pandas as pd
import h5py as h
from tqdm import tqdm

dir_path = os.getcwd()

import codebase.topo_analysis_functions as ttools
import codebase.cube.slicedice as cut

import matplotlib
from matplotlib import pyplot as plt

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))
meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')

from skimage.morphology import binary_erosion as be
from skimage.morphology import binary_dilation as bd
# from matplotlib.path import Path
# agu_data.visit(print)
def get_strat_bounary_verts(boolArr):
    sh1 = bd(boolArr) ^ boolArr
    return np.where(sh1)


chan10 = np.copy(agu_data['channelPresence10_Z'])
chan15 = np.copy(agu_data['channelPresence15_Z'])
chan30 = np.copy(agu_data['channelPresence30_Z'])

# chan10_wZero = np.copy(chan10)
# chan10_wZero[np.isnan(chan10_wZero)] = 0
# chan15_wZero = np.copy(chan15)
# chan15_wZero[np.isnan(chan15_wZero)] = 0
# chan30_wZero = np.copy(chan30)
# chan30_wZero[np.isnan(chan30_wZero)] = 0

chan10 = ma.masked_array(chan10, mask = np.isnan(chan10))
chan15 = ma.masked_array(chan15, mask = np.isnan(chan15))
chan30 = ma.masked_array(chan30, mask = np.isnan(chan30))

apex12 = (79, 87)
apex19 = (78, 78)

rs = np.array(
    [[0.5, 0.75, 1],
    [0.5, 0.9, 1.5],
    [0.5, 0.8, 1.2]]
)

slc10 = []
slc15 = []
slc30 = []

strat10 = np.copy(agu_data['synstratMinqv10'])
strat15 = np.copy(agu_data['synstratMinqv15'])
strat30 = np.copy(agu_data['synstratMinqv30'])

strat30[:, 5:11, 147:160] = np.nan

for r in rs.T:
    slc10.append(cut.circular_slice(chan10, apex12, int(200 * r[0])))
    slc15.append(cut.circular_slice(chan15, apex19, int(200 * r[1])))
    slc30.append(cut.circular_slice(chan30, apex19, int(200 * r[2])))

slcs = [slc10, slc15, slc30]

maps = [chan10.mean(axis = 0), chan15.mean(axis = 0), chan30.mean(axis = 0)]
f = plt.figure(figsize = (6, 6))

gg = f.add_gridspec(ncols = 3, nrows = 4,
    height_ratios = [3, 1, 1, 1], width_ratios = [1,1,1],
    wspace=0.15, hspace=-0.10#, top=0.95, bottom=0.05, left=0.17, right=0.845
)

xsecAsp = 8

map_ax = []

maps_saved = []

for i, top in enumerate(maps):
    if i == 0:
        map_ax.append(f.add_subplot(gg[0, i]))
    else:
        map_ax.append(f.add_subplot(gg[0, i]))
        # map_ax.append(f.add_subplot(gg[0, i], sharey = map_ax[i-1], sharex = map_ax[i-1]))
    topMin = np.nanmin(top)
    topMax = np.nanmax(top)
    maps_saved.append(map_ax[i].imshow(top,
        extent = [0, (top.shape[1] / 200), 0, (top.shape[0] / 200)],
        origin = 'lower', cmap = 'Greys', vmin = topMin, vmax = topMax
    ))

titles = ['(a) No Flooding\n($Q_v = 1$)', '(b) Low-Intensity\n($Q_v = 1.5$)', '(c) High-Intensity\n($Q_v = 3$)']

base_txtsize = 8

k = 0
for i, tit, mskdat in zip(map_ax, titles, maps):
    i.set_title(tit, fontsize = base_txtsize)
    i.spines['top'].set_visible(False)
    i.spines['bottom'].set_visible(False)
    i.spines['right'].set_visible(False)
    i.spines['left'].set_visible(False)
    i.set_xlabel('Distance (m)', fontsize = base_txtsize - 1)
    i.set_yticks([])
    i.set_xticks(np.arange(0, 2.1, 0.5))
    i.set_xticklabels(i.get_xticks(), fontsize = base_txtsize - 2)
    # i.set_yticks(np.arange(0, 2.6, 0.5))
    # i.set_yticklabels(i.get_yticks(), fontsize = base_txtsize - 2)
    for j in range(len(rs)):
        yS, xS = get_strat_bounary_verts(mskdat.mask)
        x = slcs[k][j][2]/200
        y = slcs[k][j][1]/200
        sec = slcs[k][j][0]
        nn = np.any(~np.isnan(sec), axis = 0)
        i.plot(xS/200, yS/200, 'k,', markersize=0.25)
        i.plot(x[nn], y[nn], 'r')
    k+=1

map_ax[0].set_yticks(np.arange(0, 2.6, 0.5))
map_ax[0].set_yticklabels(map_ax[0].get_yticks(), fontsize = base_txtsize - 2)
map_ax[0].set_ylabel('Distance (m)', fontsize = base_txtsize - 1)

xsec10_ax = []

secXMax = 2.5
secYMax = 0.1

for i in range(len(rs)):
#     if i == 0:
    xsec10_ax.append(f.add_subplot(gg[i+1, 0], xlim = (0, secXMax), ylim = (0, secYMax)))
    x = slcs[0][i][2]/200
    y = slcs[0][i][1]/200
    sec = slcs[0][i][0]
    nn = np.any(~np.isnan(sec), axis = 0)
    dx = np.diff(x[nn], prepend = x[nn][0])
    dy = np.diff(y[nn], prepend = y[nn][0])
    d = np.sqrt(dx * dx + dy * dy).cumsum()
    yb, xb = get_strat_bounary_verts(np.isnan(sec))
    xsec10_ax[i].imshow(sec, origin = 'lower', cmap = 'Greys',
        extent = [0, np.max(d), 0, (sec.shape[0] / 1000)]
    )
    xsec10_ax[i].spines['top'].set_visible(False)
    xsec10_ax[i].set_facecolor('#d1d1d1')
    xsec10_ax[i].set(aspect = xsecAsp)
    xsec10_ax[i].spines['bottom'].set_visible(False)
    xsec10_ax[i].spines['right'].set_visible(False)
    xsec10_ax[i].spines['left'].set_visible(False)
    xsec10_ax[i].set_yticks([0.0, 0.05])
    xsec10_ax[i].set_yticklabels(xsec10_ax[i].get_yticks()*1000, fontsize = base_txtsize - 2)
    xsec10_ax[i].set_xticks([])
xsec10_ax[-1].set_xticks(np.arange(0, 2.6, 0.5))
xsec10_ax[-1].set_xticklabels(xsec10_ax[i].get_xticks(), fontsize = base_txtsize - 2)
xsec10_ax[-1].set_xlabel('Distance (m)', fontsize = base_txtsize - 1)
xsec10_ax[1].set_ylabel('Elevation (m)', fontsize = base_txtsize - 1)

xsec15_ax = []

for i in range(len(rs)):
    xsec15_ax.append(f.add_subplot(gg[i+1, 1], xlim = (0, secXMax), ylim = (0, secYMax)))
    x = slcs[1][i][2]/200
    y = slcs[1][i][1]/200
    sec = slcs[1][i][0]
    nn = np.any(~np.isnan(sec), axis = 0)
    dx = np.diff(x[nn], prepend = x[nn][0])
    dy = np.diff(y[nn], prepend = y[nn][0])
    d = np.sqrt(dx * dx + dy * dy).cumsum()
    yb, xb = get_strat_bounary_verts(np.isnan(sec))
    xsec15_ax[i].imshow(sec, origin = 'lower', cmap = 'Greys',
        extent = [0, np.max(d), 0, (sec.shape[0] / 1000)]
    )
    xsec15_ax[i].set(aspect = xsecAsp)
    xsec15_ax[i].spines['top'].set_visible(False)
    xsec15_ax[i].set_facecolor('#d1d1d1')
    xsec15_ax[i].spines['bottom'].set_visible(False)
    xsec15_ax[i].spines['right'].set_visible(False)
    xsec15_ax[i].spines['left'].set_visible(False)
    xsec15_ax[i].set_yticks([])
    xsec15_ax[i].set_xticks([])
xsec15_ax[-1].set_xticks(np.arange(0, 2.6, 0.5))
xsec15_ax[-1].set_xticklabels(xsec15_ax[i].get_xticks(), fontsize = base_txtsize - 2)
xsec15_ax[-1].set_xlabel('Distance (m)', fontsize = base_txtsize - 1)

xsec15_ax[-1].set_xlabel('Distance (m)')

xsec30_ax = []

for i in range(len(rs)):
    xsec30_ax.append(f.add_subplot(gg[i+1, 2], xlim = (0, secXMax), ylim = (0, secYMax)))
    x = slcs[2][i][2]/200
    y = slcs[2][i][1]/200
    sec = slcs[2][i][0]
    nn = np.any(~np.isnan(sec), axis = 0)
    dx = np.diff(x[nn], prepend = x[nn][0])
    dy = np.diff(y[nn], prepend = y[nn][0])
    d = np.sqrt(dx * dx + dy * dy).cumsum()
    yb, xb = get_strat_bounary_verts(np.isnan(sec))
    xsec30_ax[i].imshow(sec, origin = 'lower', cmap = 'Greys',
        extent = [0, np.max(d), 0, (sec.shape[0] / 1000)]
    )
    # xsec30_ax[i].plot(xb, yb, 'k,', markersize = 0.5)
    xsec30_ax[i].spines['top'].set_visible(False)
    xsec30_ax[i].set_facecolor('#d1d1d1')
    xsec30_ax[i].set(aspect = xsecAsp)
    xsec30_ax[i].spines['bottom'].set_visible(False)
    xsec30_ax[i].spines['right'].set_visible(False)
    xsec30_ax[i].spines['left'].set_visible(False)
    xsec30_ax[i].set_yticks([])
    xsec30_ax[i].set_xticks([])
xsec30_ax[-1].set_xticks(np.arange(0, 2.6, 0.5))
xsec30_ax[-1].set_xticklabels(xsec30_ax[i].get_xticks(), fontsize = base_txtsize - 2)
xsec30_ax[-1].set_xlabel('Distance (m)', fontsize = base_txtsize - 1)
xsec30_ax[-1].set_xlabel('Distance (m)')

xsec15_ax[1].annotate('channel', xy=(0.82, 0.035),  xycoords='data',
            xytext=(0.9, 1.3), textcoords='axes fraction',
            fontsize = 8, color = 'red',
            arrowprops=dict(
                color='r', alpha = 0.75,
                arrowstyle="-",  connectionstyle="arc3,rad=.4"
            ),
            horizontalalignment='right', verticalalignment='top',
            )

xsec10_ax[0].annotate('floodplain', xy=(1.15, 0.035),  xycoords='data',
            xytext=(0.9, -0.05), textcoords='axes fraction',
            fontsize = 8, color = 'red',
            arrowprops=dict(
                color='r', alpha = 0.75,
                arrowstyle="-",  connectionstyle="arc3,rad=-.4"
            ),
            horizontalalignment='right', verticalalignment='top',
            )

map_ax[0].text(-0.02, 0.21, 'A', fontsize = 8, color = 'red', transform=map_ax[0].transAxes)
map_ax[0].text(0.28, -0.065, 'A$^\prime$', fontsize = 8, color = 'red', transform=map_ax[0].transAxes)
map_ax[0].text(-0.065, 0.35, 'B', fontsize = 8, color = 'red', transform=map_ax[0].transAxes,
               bbox = dict(boxstyle = 'square, pad=0.1', fc = 'w', ec = 'w'))
map_ax[0].text(0.41, -0.065, 'B$^\prime$', fontsize = 8, color = 'red', transform=map_ax[0].transAxes)
map_ax[0].text(-0.07, 0.46, 'C', fontsize = 8, color = 'red', transform=map_ax[0].transAxes)
map_ax[0].text(0.51, -0.065, 'C$^\prime$', fontsize = 8, color = 'red', transform=map_ax[0].transAxes)

map_ax[1].text(-0.02, 0.21, 'D', fontsize = 8, color = 'red', transform=map_ax[1].transAxes)
map_ax[1].text(0.28, -0.065, 'D$^\prime$', fontsize = 8, color = 'red', transform=map_ax[1].transAxes)
map_ax[1].text(-0.07, 0.42, 'E', fontsize = 8, color = 'red', transform=map_ax[1].transAxes)
map_ax[1].text(0.5, 0.04, 'E$^\prime$', fontsize = 8, color = 'red', transform=map_ax[1].transAxes)
map_ax[1].text(0.05, 0.7, 'F', fontsize = 8, color = 'red', transform=map_ax[1].transAxes)
map_ax[1].text(0.67, 0.3, 'F$^\prime$', fontsize = 8, color = 'red', transform=map_ax[1].transAxes)

map_ax[2].text(-0.04, 0.21, 'G', fontsize = 8, color = 'red', transform=map_ax[2].transAxes)
map_ax[2].text(0.28, -0.065, 'G$^\prime$', fontsize = 8, color = 'red', transform=map_ax[2].transAxes)
map_ax[2].text(-0.065, 0.38, 'H', fontsize = 8, color = 'red', transform=map_ax[2].transAxes)
map_ax[2].text(0.43, -0.065, 'H$^\prime$', fontsize = 8, color = 'red', transform=map_ax[2].transAxes)
map_ax[2].text(-0.065, 0.57, 'I', fontsize = 8, color = 'red', transform=map_ax[2].transAxes)
map_ax[2].text(0.62, -0.04, 'I$^\prime$', fontsize = 8, color = 'red', transform=map_ax[2].transAxes)

cax = f.add_axes([0.92, 0.55, 0.01, 0.25])
cbar = f.colorbar(
    maps_saved[-1], cax = cax,
    ticks = np.arange(0, 1.1, 0.25)
)
cbar.ax.set_yticklabels(cbar.ax.get_yticks(), fontsize = base_txtsize-2)
cbar.ax.set_ylabel('% time spent as channel', fontsize = base_txtsize-1)

xsec10_ax[0].text(0, 0.07, 'A', fontsize = 8, color = 'red')
xsec10_ax[0].text(1.4, 0.065, 'A$^\prime$', fontsize = 8, color = 'red')
xsec10_ax[1].text(0, 0.06, 'B', fontsize = 8, color = 'red')
xsec10_ax[1].text(2.0, 0.06, 'B$^\prime$', fontsize = 8, color = 'red')
xsec10_ax[2].text(0, 0.06, 'C', fontsize = 8, color = 'red')
xsec10_ax[2].text(2.2, 0.06, 'C$^\prime$', fontsize = 8, color = 'red')

xsec15_ax[0].text(0, 0.05, 'D', fontsize = 8, color = 'red')
xsec15_ax[0].text(1.5, 0.04, 'D$^\prime$', fontsize = 8, color = 'red')
xsec15_ax[1].text(0, 0.055, 'E', fontsize = 8, color = 'red')
xsec15_ax[1].text(1.6, 0.05, 'E$^\prime$', fontsize = 8, color = 'red')
xsec15_ax[2].text(0.25, 0.035, 'F', fontsize = 8, color = 'red')
xsec15_ax[2].text(1.15, 0.035, 'F$^\prime$', fontsize = 8, color = 'red')

xsec30_ax[0].text(0, 0.07, 'G', fontsize = 8, color = 'red')
xsec30_ax[0].text(1.4, 0.065, 'G$^\prime$', fontsize = 8, color = 'red')
xsec30_ax[1].text(0, 0.055, 'H', fontsize = 8, color = 'red')
xsec30_ax[1].text(2.0, 0.055, 'H$^\prime$', fontsize = 8, color = 'red')
xsec30_ax[2].text(0, 0.05, 'I', fontsize = 8, color = 'red')
xsec30_ax[2].text(2.3, 0.05, 'I$^\prime$', fontsize = 8, color = 'red')

f.subplots_adjust(left = 0.077, bottom = 0.075, top = 0.95, right = 0.95)

f.savefig(fname = os.path.join(dir_path, 'analysis/chapter_2/figures/output/channel_abundance.pdf'),
          transparent = False, dpi = 200)

agu_data.close()
