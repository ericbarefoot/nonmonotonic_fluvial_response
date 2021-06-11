import os
import numpy as np
import pandas as pd
import h5py as h
from tqdm import tqdm
import seaborn as sns

import codebase.topo_analysis_functions as ttools
import codebase.cube.slicedice as cut

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

dir_path = os.getcwd()

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))
meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')
# agu_data.visit(print)

# strat10 = np.copy(agu_data['synstratMinqv10'])
strat15 = np.copy(agu_data['synstratMinqv15'])
strat30 = np.copy(agu_data['synstratMinqv30'])
topo15 = np.copy(agu_data['topoqv15'])
topo30 = np.copy(agu_data['topoqv30'])

f = plt.figure(figsize = (6,7), constrained_layout = True)

gg = f.add_gridspec(ncols = 2, nrows = 3, width_ratios = [1, 1], height_ratios = [2, 1, 1])

ddd = [strat15, strat30]

topo_for_pictures = [topo15, topo30]

iinst = [25, 51]

pt115 = (285, 200)
pt215 = (235, 240)

pt130 = (235, 240)
pt230 = (155, 300)

pts = [(pt115, pt215), (pt130, pt230)]

mapax = []
secax = []

titles = ['(a) Low-Intensity\n($Q_v = 1.5$)', '(b) High-Intensity\n($Q_v = 3$)']

base_fontsize = 8

xsecann = [('A', 'A$^\prime$'), ('B', 'B$^\prime$')]

telev = 0.075

lag = 25

for j, (p, ii, strat, tit, xsa, topo) in enumerate(zip(pts, iinst, ddd, titles, xsecann, topo_for_pictures)):
    mapax.append(f.add_subplot(gg[0, j]))
    m = mapax[j]
    imgtopo = topo[ii] - np.nanmin(topo[ii])
    # print(np.nanmax(topo[ii]))
    imgmap = m.imshow(
        imgtopo, cmap = 'Greys_r',
        extent = [0, (imgtopo.shape[1] / 200), 0, (imgtopo.shape[0] / 200)],
        origin = 'lower', vmin = 0, vmax = telev
    )
    m.set_title(tit, fontsize = base_fontsize, pad = -10, y = 1)
    m.spines['top'].set_visible(False)
    m.spines['bottom'].set_visible(False)
    m.spines['right'].set_visible(False)
    m.spines['left'].set_visible(False)
    m.set_xticks(np.arange(0, 2.6, 0.5))
    m.set_xticklabels(m.get_xticks(), fontsize = base_fontsize - 2)
    m.set_xlabel('Distance (m)', fontsize = base_fontsize - 1)
    m.set_yticks([])

    m.plot(
        (p[0][1] / 200, p[1][1] / 200),
        (p[0][0] / 200, p[1][0] / 200),
        'ro-', markersize = 3
    )
    # print(p)
    m.text(
        (p[0][1] / 200)-0.05, (p[0][0] / 200)+0.05, xsa[0],
        fontsize = base_fontsize - 2, color = 'w',
    )
    m.text(
        (p[1][1] / 200)+0.05, (p[1][0] / 200)-0.1, xsa[1],
        fontsize = base_fontsize - 2, color = 'w',
    )

    if j == 0:
        secax.append(f.add_subplot(gg[(j + 1), :]))
    else:
        secax.append(f.add_subplot(gg[(j + 1), :], sharey = secax[0]))

    ss = secax[j]

    slc = cut.chord_slice(strat, p[0], p[1])

    x = slc[2]/200
    y = slc[1]/200
    sec = slc[0]
    nn = np.any(~np.isnan(sec), axis = 0)
    dx = np.diff(x[nn], prepend = x[nn][0])
    dy = np.diff(y[nn], prepend = y[nn][0])
    d = np.sqrt(dx * dx + dy * dy).cumsum()

    hlsec = sec[ii:(ii+lag), nn]
    secMin = np.nanmin(hlsec)
    hlsec = hlsec - secMin

    if j == 0:
        d = d + 0.1
        pass
    else:
        hlsec = hlsec + 0.015
        pass

    segments_r = LineCollection(
        [np.column_stack([d, y]) for y in hlsec],
        alpha = 1, linewidth = 1, cmap = 'Greys', norm=plt.Normalize(vmin=-5, vmax=30)
    )
    segments_r.set_array(np.arange(5,lag+5))

    ss.add_collection(segments_r)
    ss.set_yticks([])
    ss.set_xlim((0, 0.5))
    ss.spines['top'].set_visible(False)
    ss.spines['bottom'].set_visible(False)
    ss.spines['right'].set_visible(False)
    ss.spines['left'].set_visible(False)
    ss.set_xticks([])
    ss.set_yticks(np.arange(0, 0.045, 0.005))
    ss.set_yticklabels(secax[0].get_yticks()*1000, fontsize = base_fontsize - 2)
    ss.set_ylabel('Elevation (mm)', fontsize = base_fontsize - 1)
    # ss.set_xticklabels(np.round(ss.get_xticks()*100), fontsize = base_fontsize - 2)
    # ss.set_xlabel('Distance (cm)', fontsize = base_fontsize - 1)

#     ss.set_ylim((np.min(hlsec), np.max(hlsec)))
secax[0].set_ylim((0, 0.045))
secax[1].set_ylim((0, 0.045))

secax[1].set_xticks(np.arange(0, 0.51, 0.1))
secax[1].set_xticklabels(np.round(secax[1].get_xticks()*100), fontsize = base_fontsize - 2)
secax[1].set_xlabel('Distance (cm)', fontsize = base_fontsize - 1)

cax = f.add_axes([0.92, 0.6, 0.01, 0.25])
cbar = f.colorbar(
    imgmap, cax = cax,
    ticks = np.arange(0, telev+0.001, 0.01)
)
cbar.ax.set_yticklabels(cbar.ax.get_yticks()*1000, fontsize = base_fontsize-2)
cbar.ax.set_ylabel('Elevation (mm)', fontsize = base_fontsize-1)

secax[0].text(
    0.10, 0.030, 'A', fontsize = base_fontsize,
    color = 'red', fontweight = 'bold'
)
secax[0].text(
    0.42, 0.035, 'A$^\prime$', fontsize = base_fontsize,
    color = 'red', fontweight = 'bold'
)

secax[1].text(
    0.01, 0.025, 'B', fontsize = base_fontsize,
    color = 'red', fontweight = 'bold'
)
secax[1].text(
    0.49, 0.035, 'B$^\prime$', fontsize = base_fontsize,
    color = 'red', fontweight = 'bold'
)

secax[0].annotate('', xy=(0.25, 0.02),  xycoords='data',
            xytext=(0.25, 0.005), textcoords='data',
            fontsize = base_fontsize - 2, color = 'r',
            arrowprops=dict(
                color='r', alpha = 1,
                arrowstyle="->",  connectionstyle="arc3,rad=0",
                relpos = (1, 0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

secax[0].annotate('levee aggradation', xy=(0.2, 0.034),  xycoords='data',
            xytext=(0.25, 0.045), textcoords='data',
            fontsize = base_fontsize - 1, color = 'r',
            arrowprops=dict(
                color='r', alpha = 1,
                arrowstyle="->",  connectionstyle="arc3,rad=0.2",
                relpos = (0, 0)
            ),
            horizontalalignment='center', verticalalignment='top',
            )
secax[0].annotate('levee aggradation', xy=(0.30, 0.033),  xycoords='data',
            xytext=(0.25, 0.045), textcoords='data',
            fontsize = base_fontsize - 1, color = 'r',
            arrowprops=dict(
                color='r', alpha = 1,
                arrowstyle="->",  connectionstyle="arc3,rad=-0.2",
                relpos = (1, 0)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

secax[0].annotate('', xy=(0.21, 0.025),  xycoords='data',
            xytext=(0.21, 0.00), textcoords='data',
            fontsize = base_fontsize - 2, color = 'r',
            arrowprops=dict(
                color='r', alpha = 1,
                arrowstyle="<->",
                connectionstyle="arc,angleA=180,angleB=180,armA=15,armB=15,rad=2",
                relpos = (1, 0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

secax[0].text(
    0.17, 0.0125, 'deep channel', fontsize = base_fontsize,
    color = 'red', rotation = 'vertical', va = 'center'
)

secax[1].annotate('', xy=(0.45, 0.033),  xycoords='data',
            xytext=(0.23, 0.029), textcoords='data',
            fontsize = base_fontsize - 2, color = 'r',
            arrowprops=dict(
                color='r', alpha = 1,
                arrowstyle="->",  connectionstyle="arc3,rad=-0.4",
                relpos = (1, 0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

secax[1].annotate('', xy=(0.3, 0.017),  xycoords='data',
            xytext=(0.14, 0.017), textcoords='data',
            fontsize = base_fontsize - 2, color = 'r',
            arrowprops=dict(
                color='r', alpha = 1,
                arrowstyle="<->",
                connectionstyle="arc,angleA=-90,angleB=-90,armA=15,armB=15,rad=2",
                relpos = (1, 0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )
secax[1].text(
    0.22, 0.007, 'wide, shallow channel', fontsize = base_fontsize,
    color = 'red', rotation = 'horizontal', va = 'center', ha = 'center'
)

secax[1].text(
    0.25, 0.038, 'avulsion', fontsize = base_fontsize,
    color = 'red', rotation = 35, va = 'center', ha = 'center'
)

secax[1].annotate('', xy=(0.21, 0.026),  xycoords='data',
            xytext=(0.16, 0.020), textcoords='data',
            fontsize = base_fontsize - 2, color = 'r',
            arrowprops=dict(
                color='r', alpha = 1,
                arrowstyle="->",  connectionstyle="arc3,rad=0",
                relpos = (1, 0.5)
            ),
            horizontalalignment='center', verticalalignment='top',
            )

secax[0].annotate(
    'pulses of in-channel aggradation\nbalanced by levee growth.',
    xytext=(0.32, 0.013), textcoords='data',
    xy=(0.265, 0.013), xycoords='data',
    fontsize = base_fontsize, color = 'r',
    arrowprops=dict(
        color='r', alpha = 1,
        arrowstyle="->",
        connectionstyle="arc3,rad=0",
        relpos = (0, 0.5)
    ),
    horizontalalignment='left', verticalalignment='center',
)

secax[1].annotate(
    'in-channel aggradation\nand lateral accretion\ndestabilizes channel',
    xytext=(0.14, 0.035), textcoords='data',
    xy=(0.18, 0.0275), xycoords='data',
    fontsize = base_fontsize, color = 'r',
    arrowprops=dict(
        color='r', alpha = 1,
        arrowstyle="->",
        connectionstyle="arc3,rad=-0.2",
        relpos = (1, 0.5)
    ),
    horizontalalignment='right', verticalalignment='center',
)

# 0.28, 0.035,
# arc,angleA=0,angleB=-25,armA=30,armB=70,rad=5

mapax[0].set_yticks(np.arange(0, 2.6, 0.5))
mapax[0].set_yticklabels(mapax[0].get_yticks(), fontsize = base_fontsize - 2)
mapax[0].set_ylabel('Distance (m)', fontsize = base_fontsize - 1)

f.savefig(fname = os.path.join(dir_path, 'analysis/chapter_2/figures/output/channel_evolution.pdf'),
          transparent = False, dpi = 200)
