
import os
import numpy as np
import pandas as pd
import h5py as h
from tqdm import tqdm
import matplotlib
from matplotlib import pyplot as plt

dir_path = os.getcwd()

import sys
sys.path.insert(0, dir_path)

import codebase.cube.slicedice as cut

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')

meta = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))

timepick = 54
tind = agu_data['topoqv30'].attrs['IDs'][timepick]
iind = int(np.where(agu_data['imagRefqv30'].attrs['IDs'] == (tind - 4))[0])
mind = int(np.where(agu_data['channelMasksqv30'].attrs['IDs'] == (tind - 4))[0])
mmind = int(np.where(agu_data['manualChannelMasksqv30'].attrs['IDs'] == (tind - 4))[0])

apex19 = (78, 78)

qv = 1.5

r = int(1.2 * 200)

mstr = 'manualChannelMasksqv30'

mask = np.copy(agu_data['manualChannelMasksqv30'])
imag = np.copy(agu_data['imagRefqv30'])
topo = np.copy(agu_data['topoqv30'])

sec, x, y = cut.circular_slice(topo, apex19, r)
nn = np.any(~np.isnan(sec), axis = 0)
dx = np.diff(x[nn], prepend = x[nn][0])
dy = np.diff(y[nn], prepend = y[nn][0])
d = np.sqrt(dx * dx + dy * dy).cumsum()
maskI = mask   #[match, :]
mskSlc = maskI[:, x[nn], y[nn]]
mskSlc[:, -1] = False
mskDiff = np.diff(mskSlc.astype('int'), axis = -1, prepend = 0)
rbs = np.where(mskDiff == 1)[1] - 10
rbsT = np.where(mskDiff == 1)[0]
lbs = np.where(mskDiff == -1)[1] + 0
lbsT = np.where(mskDiff == -1)[0]

matching_time = rbsT == timepick

#################################
# fiftyfour = rbsT == 54
#
# # plt.clf()
# # plt.imshow(mskDiff)
# # plt.plot(rbs[fiftyfour], rbsT[fiftyfour], 'r.')
# # plt.plot(lbs[fiftyfour], lbsT[fiftyfour], 'r.')
#
# rb, lb = (rbs[fiftyfour], lbs[fiftyfour])
#
# plt.clf()
# plt.imshow(mask[mmind])
# plt.plot(y[nn], x[nn], 'r')
# plt.plot(y[lb], x[lb], 'w.')
# plt.plot(y[rb], x[rb], 'w.')
#

#################################

chanDimData = {k: [] for k in ['qv', 'time', 'radius', 'width', 'depth', 'xrb', 'yrb', 'xlb', 'ylb', 'tID']}

for rb, lb, t in zip(rbs[matching_time], lbs[matching_time], rbsT[matching_time]):
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
chanpts = chandat[chandat.tID == (tind - 4)]

f = plt.figure(figsize = (6,4))

gg = f.add_gridspec(
    ncols = 3, nrows = 2,
    width_ratios = [1,1,1], height_ratios = [2,1],
    wspace=0.15, hspace=0.5
)

dd = [imag[iind], mask[mmind], topo[timepick]]
tt = ['(a) Image', '(b) Channel Mask', '(c) Topography']

img = f.add_subplot(gg[0, 0])
msk = f.add_subplot(gg[0, 1])
man = f.add_subplot(gg[0, 2])
xsec = f.add_subplot(gg[1, :])

ax = (img, msk, man)

base_fontsize = 8

img.set_ylabel('Distance (m)', fontsize = base_fontsize - 2)

j = 0
for i, dif, tit in zip(ax[0:3], dd, tt):
    # i.set_facecolor('#000000')
    i.set_title(tit, fontsize = base_fontsize)
    i.spines['top'].set_visible(False)
    i.spines['bottom'].set_visible(False)
    i.spines['right'].set_visible(False)
    i.spines['left'].set_visible(False)
    i.plot(y[nn] / 200, x[nn] / 200, 'r-')
    if j == 0:
        i.imshow(dif,
            extent = [0, (dif.shape[1] / 200), 0, (dif.shape[0] / 200)],
            origin = 'lower')
    else:
        imgmap = i.imshow(dif,
            extent = [0, (dif.shape[1] / 200), 0, (dif.shape[0] / 200)],
            origin = 'lower', cmap = 'Greys_r')
    j += 1
    i.plot(chanpts.yrb / 200, chanpts.xrb / 200, 'w.', markersize = 2)
    i.plot(chanpts.ylb / 200, chanpts.xlb / 200, 'w.', markersize = 2)
    i.set_xticks(np.arange(0, 2.6, 0.5))
    i.set_xticklabels(i.get_xticks(), fontsize = base_fontsize - 2)
    i.set_xlabel('Distance (m)', fontsize = base_fontsize - 1)
    i.set_yticks([])

ax[2].text(1.35, 1.2, 'A', fontsize = base_fontsize - 2, color = 'w')
ax[2].text(1.6, 0.75, 'A$^\prime$', fontsize = base_fontsize - 2, color = 'w')

ax[0].set_yticks(np.arange(0, 2.6, 0.5))
ax[0].set_yticklabels(ax[0].get_yticks(), fontsize = base_fontsize - 2)
ax[0].set_ylabel('Distance (m)', fontsize = base_fontsize - 1)

for i in range(chanpts.shape[0]):
    slc, xx, yy = cut.chord_slice(
        topo,
        (int(chanpts.iloc[i].xrb), int(chanpts.iloc[i].yrb)),
        (int(chanpts.iloc[i].xlb), int(chanpts.iloc[i].ylb))
    )
    dx = np.diff(xx, prepend = xx[0])
    dy = np.diff(yy, prepend = yy[0])
    d = np.sqrt(dx * dx + dy * dy).cumsum()
    xsec.plot(d / 200, slc[timepick], 'k')
xsec.set_title('(d) Channel Cross-section', fontsize = base_fontsize)
xsec.spines['top'].set_visible(False)
xsec.spines['bottom'].set_visible(False)
xsec.spines['right'].set_visible(False)
xsec.spines['left'].set_visible(False)
xsec.set_xlabel('Distance (m)', fontsize = base_fontsize - 2)
xsec.set(aspect = 15)
xsec.yaxis.tick_right()
xsec.yaxis.set_label_position("right")
xsec.set_xticks(np.arange(0, 0.21, 0.05))
xsec.set_xticklabels(np.round(xsec.get_xticks(), 2)*100, fontsize = base_fontsize - 2)
xsec.set_xlabel('Distance (cm)', fontsize = base_fontsize - 1)
xsec.set_yticks(np.arange(0.08, 0.087, 0.002))
xsec.set_yticklabels(np.round(xsec.get_yticks(), 3) * 1000, fontsize = base_fontsize - 2)
xsec.set_ylabel('Elevation (mm)', fontsize = base_fontsize - 1)

xsec.text(-0.005, 0.0845, 'A', fontsize = base_fontsize, color = 'red', fontweight = 'bold')
xsec.text(0.21, 0.085, 'A$^\prime$', fontsize = base_fontsize, color = 'red', fontweight = 'bold')


f.subplots_adjust(left = 0.075)


ax[2].annotate('', xy=(1.6, 0.75),  xycoords='data',
            xytext=(0.75, 0.2), textcoords='figure fraction',
            fontsize = 8, color = 'r',
            arrowprops=dict(
                color='r', alpha = 0.75,
                arrowstyle="<-",
                connectionstyle="arc,angleA=0,angleB=-25,armA=30,armB=70,rad=5",
                relpos = (1, 0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

cax = f.add_axes([0.87, 0.6, 0.01, 0.25])
cbar = f.colorbar(
    imgmap, cax = cax,
    ticks = np.arange(0.070, 0.161, 0.02)
)
cbar.ax.set_yticklabels((cbar.ax.get_yticks()*1000).round(), fontsize = base_fontsize-2)
cbar.ax.set_ylabel('Elevation (mm)', fontsize = base_fontsize-1)

f.savefig(fname = os.path.join(dir_path, 'analysis/chapter_2/figures/output/channel_geom_procedure.pdf'), transparent = False, dpi = 200)


# xsec30_ax[-1].set_xticks(np.arange(0, 2.6, 0.5))
# xsec30_ax[-1].set_xticklabels(xsec15_ax[i].get_xticks(), fontsize = base_txtsize - 2)
# xsec30_ax[-1].set_xlabel('Distance (m)', fontsize = base_txtsize - 1)

agu_data.close()
