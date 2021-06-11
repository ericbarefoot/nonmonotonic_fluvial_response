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
from matplotlib.collections import LineCollection

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))
meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')

from skimage.morphology import binary_erosion as be
from matplotlib.path import Path

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

lwd = 0.5

backWallROI12 = [
    (159.65764557773477, 2.512843820692325),
    (160.9917367528929, 8.38284499138814),
    (8.63852454983325, 162.07014836960593),
    (2.234886909074177, 162.33696660463755),
    (2.501705144105806, 155.1328742587836),
    (151.91991676181755, 3.313298525787218)
]
backWallROI19 = [
    (147.42873693244786, 2.1601899486271634),
    (159.97708031372878, 2.4390420237667456),
    (161.37134068942663, 8.852639751976994),
    (8.281551437799365, 162.77898522902296),
    (2.1468057847286843, 157.48079580137102),
    (4.098770310705717, 146.32671279578796)
]

topo10 = np.copy(agu_data['topoqv10'])
topo15 = np.copy(agu_data['topoqv15'])
topo30 = np.copy(agu_data['topoqv30'])

topo30[:, 5:11, 147:160] = np.nan

diff10 = topo10[55] - topo10[20]
diff15 = topo15[55] - topo15[20]
diff30 = topo30[55] - topo30[20]

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
    slc10.append(cut.circular_slice(strat10, apex12, int(200 * r[0])))
    slc15.append(cut.circular_slice(strat15, apex19, int(200 * r[1])))
    slc30.append(cut.circular_slice(strat30, apex19, int(200 * r[2])))

slcs = [slc10, slc15, slc30]

maps = [diff10, diff15, diff30]

def get_poly_forcut(sec, i1, i2, xarange):
    top = sec[[i2], nn]
    tnan = np.isnan(top)
    btm = sec[[i1], nn]
    bnan = np.isnan(btm)
    topinterp = np.interp(xarange, d[~tnan], top[~tnan])
    btminterp = np.interp(xarange, d[~bnan], btm[~bnan])
    qwe = np.hstack([np.flip(xarange), xarange])
    ewq = np.hstack([np.flip(btminterp), topinterp])
    return qwe, ewq

f = plt.figure(figsize = (6, 6))
gg = f.add_gridspec(ncols = 3, nrows = 4, height_ratios = [3, 1, 1, 1],
                    wspace=0.15, hspace=0
                   )

map_ax = []

maps_saved = []

for i, top in enumerate(maps):
    if i == 0:
        map_ax.append(f.add_subplot(gg[0, i]))
    else:
        map_ax.append(f.add_subplot(gg[0, i]))
#         map_ax.append(f.add_subplot(gg[0, i], sharey = map_ax[i-1], sharex = map_ax[i-1]))
    topMin = 0
    topMax = 20
    maps_saved.append(map_ax[i].imshow(top * 1000,
        extent = [0, (top.shape[1] / 200), 0, (top.shape[0] / 200)],
        origin = 'lower', cmap = 'Greys', vmin = topMin, vmax = topMax
    ))

titles = ['(a) No Flooding\n($Q_v = 1$)', '(b) Low-Intensity\n($Q_v = 1.5$)', '(c) High-Intensity\n($Q_v = 3$)']

base_txtsize = 8

msks = ['topsetMasksqv10', 'topsetMasksqv15', 'topsetMasksqv30']
rois = [backWallROI12, backWallROI19, backWallROI19]

k = 0
for i, tit, mskdat, roii in zip(map_ax, titles, maps, rois):
    i.set_title(tit, fontsize = base_txtsize)
    i.spines['top'].set_visible(False)
    i.spines['bottom'].set_visible(False)
    i.spines['right'].set_visible(False)
    i.spines['left'].set_visible(False)
    i.set_xlabel('Distance (m)', fontsize = base_txtsize - 1)
    i.set_xticks(np.arange(0, 2.1, 0.5))
    i.set_xticklabels(i.get_xticks(), fontsize = base_txtsize - 2)
    i.set_yticks([])
#     i.set_yticks(np.arange(0, 2.6, 0.5))
#     i.set_yticklabels(i.get_yticks(), fontsize = base_txtsize - 2)
    for j in range(len(rs)):
        yS, xS = get_shore_verts(~np.isnan(mskdat), roii, 6, 6)
        x = slcs[k][j][2]/200
        y = slcs[k][j][1]/200
        sec = slcs[k][j][0]
        nn = np.any(~np.isnan(sec), axis = 0)
        i.plot(xS/200, yS/200, 'k.', markersize=0.25)
        i.plot(x[nn], y[nn], 'r')
    k+=1
map_ax[0].set_yticks(np.arange(0, 2.6, 0.5))
map_ax[0].set_yticklabels(map_ax[0].get_yticks(), fontsize = base_txtsize - 2)
map_ax[0].set_ylabel('Distance (m)', fontsize = base_txtsize - 1)

secXMax = 2.5
secYMax = 0.075

xsec10_ax = []

for i in range(len(rs)):
    xsec10_ax.append(f.add_subplot(gg[i+1, 0], xlim = (0, secXMax), ylim = (0, secYMax)))
    x = slcs[0][i][2]/200
    y = slcs[0][i][1]/200
    sec = slcs[0][i][0]
    secMin = np.nanmin(sec)
    sec = sec - secMin
    nn = np.any(~np.isnan(sec), axis = 0)
    dx = np.diff(x[nn], prepend = x[nn][0])
    dy = np.diff(y[nn], prepend = y[nn][0])
    d = np.sqrt(dx * dx + dy * dy).cumsum()
    dMin = np.nanmin(d)
    dMax = np.nanmax(d)

    xbg, ybg = get_poly_forcut(sec, 0, -1, np.arange(dMin, dMax, 0.005))
    xfg, yfg = get_poly_forcut(sec, 20, 51, np.arange(dMin, dMax, 0.005))

    segments = LineCollection(
        [np.column_stack([d, y]) for y in sec[:, nn]],
        alpha = 0.5, color = 'k', linewidth = lwd
    )
    xsec10_ax[i].fill(xbg, ybg, fc = '#bdbdbd', ec = '#bdbdbd')
    xsec10_ax[i].fill(xfg, yfg, fc = 'r', ec = 'r')
#     xsec10_ax[i].add_collection(segments)
    xsec10_ax[i].spines['top'].set_visible(False)
    xsec10_ax[i].spines['bottom'].set_visible(False)
    xsec10_ax[i].spines['right'].set_visible(False)
    xsec10_ax[i].spines['left'].set_visible(False)
    xsec10_ax[i].set_yticks([0.0, 0.05])
    xsec10_ax[i].set_yticklabels(xsec10_ax[i].get_yticks()*1000, fontsize = base_txtsize - 2)
    xsec10_ax[i].set_xticks([])
xsec10_ax[-1].set_xticks(np.arange(0, 2.6, 0.5))
xsec10_ax[-1].set_xticklabels(xsec10_ax[i].get_xticks(), fontsize = base_txtsize - 2)
xsec10_ax[-1].set_xlabel('Distance (m)', fontsize = base_txtsize - 1)
xsec10_ax[1].set_ylabel('Elevation (mm)', fontsize = base_txtsize - 1)

xsec15_ax = []

for i in range(len(rs)):
    xsec15_ax.append(f.add_subplot(gg[i+1, 1], xlim = (0, secXMax), ylim = (0, secYMax)))
    x = slcs[1][i][2]/200
    y = slcs[1][i][1]/200
    sec = slcs[1][i][0]
    secMin = np.nanmin(sec)
    sec = sec - secMin
    nn = np.any(~np.isnan(sec), axis = 0)
    dx = np.diff(x[nn], prepend = x[nn][0])
    dy = np.diff(y[nn], prepend = y[nn][0])
    d = np.sqrt(dx * dx + dy * dy).cumsum()
    dMin = np.nanmin(d)
    dMax = np.nanmax(d)

    xbg, ybg = get_poly_forcut(sec, 0, -1, np.arange(dMin, dMax, 0.005))
    xfg, yfg = get_poly_forcut(sec, 20, 51, np.arange(dMin, dMax, 0.005))

    segments = LineCollection(
        [np.column_stack([d, y]) for y in sec[:, nn]],
        alpha = 0.5, color = 'k', linewidth = lwd
    )

    secMin = np.nanmin(sec)
    secMax = np.nanmax(sec)
    xsec15_ax[i].fill(xbg, ybg, fc = '#bdbdbd', ec = '#bdbdbd')
    xsec15_ax[i].fill(xfg, yfg, fc = 'r', ec = 'r')
#     xsec15_ax[i].add_collection(segments)
    xsec15_ax[i].set_yticks([])
    xsec15_ax[i].set_xticks([])
    xsec15_ax[i].spines['top'].set_visible(False)
    xsec15_ax[i].spines['bottom'].set_visible(False)
    xsec15_ax[i].spines['right'].set_visible(False)
    xsec15_ax[i].spines['left'].set_visible(False)
xsec15_ax[-1].set_xticks(np.arange(0, 2.6, 0.5))
xsec15_ax[-1].set_xticklabels(xsec15_ax[i].get_xticks(), fontsize = base_txtsize - 2)
xsec15_ax[-1].set_xlabel('Distance (m)', fontsize = base_txtsize - 1)

xsec30_ax = []

for i in range(len(rs)):
    xsec30_ax.append(f.add_subplot(gg[i+1, 2], xlim = (0, secXMax), ylim = (0, secYMax)))
    x = slcs[2][i][2]/200
    y = slcs[2][i][1]/200
    sec = slcs[2][i][0]
    secMin = np.nanmin(sec)
    sec = sec - secMin
    nn = np.any(~np.isnan(sec), axis = 0)
    dx = np.diff(x[nn], prepend = x[nn][0])
    dy = np.diff(y[nn], prepend = y[nn][0])
    d = np.sqrt(dx * dx + dy * dy).cumsum()
    dMin = np.nanmin(d)
    dMax = np.nanmax(d)

    xbg, ybg = get_poly_forcut(sec, 0, -1, np.arange(dMin, dMax, 0.005))
    xfg, yfg = get_poly_forcut(sec, 20, 51, np.arange(dMin, dMax, 0.005))

    segments = LineCollection(
        [np.column_stack([d, y]) for y in sec[:, nn]],
        alpha = 0.5, color = 'k', linewidth = lwd
    )
    xsec30_ax[i].fill(xbg, ybg, fc = '#bdbdbd', ec = '#bdbdbd')
    xsec30_ax[i].fill(xfg, yfg, fc = 'r', ec = 'r')

#     xsec30_ax[i].add_collection(segments)
    xsec30_ax[i].set_yticks([])
    xsec30_ax[i].set_xticks([])
    xsec30_ax[i].spines['top'].set_visible(False)
    xsec30_ax[i].spines['bottom'].set_visible(False)
    xsec30_ax[i].spines['right'].set_visible(False)
    xsec30_ax[i].spines['left'].set_visible(False)
xsec30_ax[-1].set_xticks(np.arange(0, 2.6, 0.5))
xsec30_ax[-1].set_xticklabels(xsec15_ax[i].get_xticks(), fontsize = base_txtsize - 2)
xsec30_ax[-1].set_xlabel('Distance (m)', fontsize = base_txtsize - 1)


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

cax = f.add_axes([0.92, 0.6, 0.01, 0.25])
cbar = f.colorbar(
    maps_saved[-1], cax = cax,
    ticks = np.arange(topMin, topMax+1, 5)
)
cbar.ax.set_yticklabels(cbar.ax.get_yticks(), fontsize = base_txtsize-2)
cbar.ax.set_ylabel('35hr-Accumulation (mm)', fontsize = base_txtsize-1)


xsec10_ax[0].text(-0.1, 0.037, 'A', fontsize = 8, color = 'red')
xsec10_ax[0].text(1.65, 0.037, 'A$^\prime$', fontsize = 8, color = 'red')
xsec10_ax[1].text(-0.1, 0.037, 'B', fontsize = 8, color = 'red')
xsec10_ax[1].text(2.1, 0.037, 'B$^\prime$', fontsize = 8, color = 'red')
xsec10_ax[2].text(-0.1, 0.037, 'C', fontsize = 8, color = 'red')
xsec10_ax[2].text(2.4, 0.037, 'C$^\prime$', fontsize = 8, color = 'red')

xsec15_ax[0].text(-0.1, 0.037, 'D', fontsize = 8, color = 'red')
xsec15_ax[0].text(1.6, 0.037, 'D$^\prime$', fontsize = 8, color = 'red')
xsec15_ax[1].text(-0.1, 0.037, 'E', fontsize = 8, color = 'red')
xsec15_ax[1].text(1.6, 0.037, 'E$^\prime$', fontsize = 8, color = 'red')
xsec15_ax[2].text(0.0, 0.037, 'F', fontsize = 8, color = 'red')
xsec15_ax[2].text(1.75, 0.037, 'F$^\prime$', fontsize = 8, color = 'red')

xsec30_ax[0].text(-0.1, 0.04, 'G', fontsize = 8, color = 'red')
xsec30_ax[0].text(1.6, 0.035, 'G$^\prime$', fontsize = 8, color = 'red')
xsec30_ax[1].text(-0.1, 0.037, 'H', fontsize = 8, color = 'red')
xsec30_ax[1].text(2.0, 0.037, 'H$^\prime$', fontsize = 8, color = 'red')
xsec30_ax[2].text(-0.1, 0.037, 'I', fontsize = 8, color = 'red')
xsec30_ax[2].text(2.5, 0.037, 'I$^\prime$', fontsize = 8, color = 'red')

xsec15_ax[2].annotate('delta lobe', xy=(1.4, 0.02),  xycoords='data',
            xytext=(0.3, 0.9), textcoords='axes fraction',
            fontsize = 8, color = 'k',
            arrowprops=dict(
                color='k', alpha = 0.75,
                arrowstyle="-",  connectionstyle="arc3,rad=-.4",
                relpos = (1, 0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )
xsec15_ax[1].annotate('levee', xy=(1.0, 0.03),  xycoords='data',
            xytext=(0.3, 0.9), textcoords='axes fraction',
            fontsize = 8, color = 'k',
            arrowprops=dict(
                color='k', alpha = 0.75,
                arrowstyle="-",  connectionstyle="arc3,rad=-.4",
                relpos = (1, 0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

xsec30_ax[1].annotate('dispersed\nsediment', xy=(1.6, 0.025),  xycoords='data',
            xytext=(0.92, 0.38), textcoords='figure fraction',
            fontsize = 8, alpha = 1, color = 'k',
            arrowprops=dict(
                color='k', alpha = 0.75,
                arrowstyle="-",  connectionstyle="arc3,rad=.4",
                relpos = (0,0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

xsec10_ax[1].annotate('lobe', xy=(0.3, 0.025),  xycoords='data',
            xytext=(0.3, 0.9), textcoords='axes fraction',
            fontsize = 8, color = 'k',
            arrowprops=dict(
                color='k', alpha = 0.75,
                arrowstyle="-",  connectionstyle="arc3,rad=.4",
                relpos = (0, 0.25)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

f.subplots_adjust(left = 0.077, bottom = 0.075, top = 0.95, right = 0.95)

f.savefig(fname = os.path.join(dir_path, 'analysis/chapter_2/figures/output/strat_evenness.pdf'), transparent = False, dpi = 200)

agu_data.close()
