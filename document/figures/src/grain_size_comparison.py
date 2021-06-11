
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

import numpy.ma as ma
from skimage.morphology import binary_dilation as bd

def get_strat_bounary_verts(boolArr):
    sh1 = bd(boolArr) ^ boolArr
    return np.where(sh1)

meta_data_19 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb19_metadata.csv'))
meta_data_12 = pd.read_csv(os.path.join(dir_path, 'data/raw_data/tdb12_metadata.csv'))

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')

grain_size_input_mix = pd.read_csv(
    os.path.join(dir_path, 'analysis/chapter_2/data/grainsize_curves_input_mix.csv')
)
grain_size_select = pd.read_csv(
    os.path.join(dir_path, 'analysis/chapter_2/data/grainsize_curves_examples.csv')
)
grain_size_all = pd.read_csv(
    os.path.join(dir_path, 'analysis/chapter_2/data/grainsize_curves_to_compare.csv')
)

chan10 = np.copy(agu_data['channelPresence10_Z'])
chan15 = np.copy(agu_data['channelPresence15_Z'])
chan30 = np.copy(agu_data['channelPresence30_Z'])

chan10 = ma.masked_array(chan10, mask = np.isnan(chan10))
chan15 = ma.masked_array(chan15, mask = np.isnan(chan15))
chan30 = ma.masked_array(chan30, mask = np.isnan(chan30))

apex12 = (79, 87)
apex19 = (78, 78)

# rs = np.array(
#     [
#         [0.5, 0.75, 1],
#         [0.5, 0.9, 1.5],
#         [0.5, 0.8, 1.2]
#     ]
# )

rs = np.array(
    [
        np.linspace(0.05, 1.5, 15),
        np.linspace(0.05, 1.5, 15),
        np.linspace(0.05, 1.5, 15)
    ]
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

chan_abund_data = {k: [] for k in ['r', 'p', 'qv']}

qvs = [1, 1.5, 3]

for rqv, ss, qv in zip(rs, slcs, qvs):
    for ri, r in enumerate(rqv):
        chh = np.nansum(ss[ri][0])
        tot = ss[ri][0].size
        chan_abund_data['p'].append(chh / tot)
        chan_abund_data['r'].append(r)
        chan_abund_data['qv'].append(qv)
#         print(chh / tot)
#         print(np.nansum(ss[ri][0]))
#         print(ss[ri][0].size)

ch_dat = pd.DataFrame(chan_abund_data)

ch_dat.to_csv('/home/eric/Dropbox/research/barefoot_dissertation/nonmonotonic_chapter/data/derived_data/channel_abundance_data.csv', index = False)

input_mix = grain_size_input_mix[(grain_size_input_mix['facies'] == 'channel')]
#& (grain_size_input_mix['facies'] == 'channel')

all_chan = grain_size_all[grain_size_all['facies'] == 'channel']
all_floo = grain_size_all[grain_size_all['facies'] == 'levee']

ex_chan = grain_size_select[grain_size_select['facies'] == 'channel']
ex_floo = grain_size_select[grain_size_select['facies'] == 'levee']

base_txtsize = 8

f = plt.figure(figsize = (6,3))

gg = f.add_gridspec(ncols = 4, nrows = 2,
                    height_ratios = [1, 1], width_ratios = [3,2,2,2],
                    wspace=0.4#, hspace=0.1
                   )

ch = f.add_subplot(gg[0, 0])
fp = f.add_subplot(gg[1, 0], sharey = ch)

ax =  (ch, fp)

all_chan.groupby('sample_name').plot(
    x = 'size', y = 'percent', color = 'grey', ax = ch, legend=False, linewidth = 1
)

all_floo.groupby('sample_name').plot(
    x = 'size', y = 'percent', color = 'grey', ax = fp, legend=False, linewidth = 1
)

for a in ax:
    a.plot(input_mix['size'], input_mix['avg_percent'], color = 'k', linewidth = 2)
    a.set_xscale('log')
    a.set_yticks(np.arange(0, 0.032, 0.01))
    a.set_yticklabels(a.get_yticks(), fontsize = base_txtsize - 2)
    a.set_ylabel('Percent %', fontsize = base_txtsize - 1)
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)
    pass

ex_chan.groupby('sample_name').plot(
    x = 'size', y = 'percent', color = 'red', ax = ch, legend=False, linewidth = 2
)

ex_floo.groupby('sample_name').plot(
    x = 'size', y = 'percent', color = 'red', ax = fp, legend=False, linewidth = 2
)

ch.set_xticks([])

fp.tick_params(axis = 'x', which = 'major', labelsize = base_txtsize - 2)
ch.tick_params(axis = 'x', which = 'both', length =  0)

fp.set_xlabel('Grain Size ($\mu$m)', fontsize = base_txtsize - 2)
ch.set_xlabel('', fontsize = base_txtsize - 1)

fp.set_title(
    '(b) Grain Size Distribution for\nFloodplain Deposits',
    fontsize = base_txtsize - 1, pad = -5, y = 1
)
ch.set_title(
    '(a) Grain Size Distribution for\nChannel Deposits',
    fontsize = base_txtsize - 1, pad = 0, y = 1
)

chanabund = np.copy(agu_data['channelPresence30_Z'])

chanabund = ma.masked_array(chanabund, mask = np.isnan(chanabund))

apex19 = (78, 78)

# strat = np.copy(agu_data['synstratMinqv30'])

slc = cut.circular_slice(chanabund, apex19, int(200 * 0.9))
xsecAsp = 5
secXMax = 2.5
secYMax = 0.1

xsec = f.add_subplot(gg[0, 1:], xlim = (0, secXMax), ylim = (0, secYMax))

x = slc[2]/200
y = slc[1]/200
sec = slc[0]
nn = np.any(~np.isnan(sec), axis = 0)
dx = np.diff(x[nn], prepend = x[nn][0])
dy = np.diff(y[nn], prepend = y[nn][0])
d = np.sqrt(dx * dx + dy * dy).cumsum()
yb, xb = get_strat_bounary_verts(np.isnan(sec))
xsec.imshow(sec, origin = 'lower', cmap = 'Greys',
    extent = [0, np.max(d), 0, (sec.shape[0] / 1000)]
)
xsec.set(aspect = xsecAsp)
xsec.spines['top'].set_visible(False)
xsec.set_facecolor('#d1d1d1')
xsec.spines['bottom'].set_visible(False)
xsec.spines['right'].set_visible(False)
xsec.spines['left'].set_visible(False)
xsec.set_yticks([0.0, 0.05, 0.1])
xsec.set_ylabel('Elevation (m)', fontsize = base_txtsize - 2)
xsec.set_yticklabels(xsec.get_yticks(), fontsize = base_txtsize - 2)
xsec.set_xticks(np.arange(0, 2.6, 0.5))
xsec.set_xticklabels(xsec.get_xticks(), fontsize = base_txtsize - 2)
xsec.set_xlabel('Distance (m)', fontsize = base_txtsize - 2)

chabund_ax = []

chabund_ax.append(f.add_subplot(gg[1, 1], xlim = (0, secXMax), ylim = (0, secYMax)))
chabund_ax.append(f.add_subplot(gg[1, 2], xlim = (0, secXMax), ylim = (0, secYMax)))
chabund_ax.append(f.add_subplot(gg[1, 3], xlim = (0, secXMax), ylim = (0, secYMax)))

titles = ['(d) No Flooding\n($Q_v = 1$)', '(e) Low-Intensity\n($Q_v = 1.5$)', '(f) High-Intensity\n($Q_v = 3$)']

for a, qv, ttl in zip(chabund_ax, qvs, titles):
    mtch = ch_dat['qv'] == qv
    # print(ch_dat['p'][mtch].to_numpy())
    a.plot(
        ch_dat['r'][mtch].to_numpy(),
        ch_dat['p'][mtch].to_numpy(),
        '-ok', markersize = 3
    )
    a.spines['top'].set_visible(False)
#     a.spines['bottom'].set_visible(False)
    a.set_facecolor('#00000000')
    a.spines['left'].set_visible(False)
    a.spines['right'].set_visible(False)
    a.set_title(ttl, fontsize = base_txtsize - 1, pad = -10, y = 1)
    a.set_yticks([])
    a.set_xticks(np.arange(0, 2.1, 0.5))
    a.set_xticklabels(a.get_xticks(), fontsize = base_txtsize - 2)
    a.set_xlabel('Radial Distance (m)', fontsize = base_txtsize - 2)
    a.set_xlim((0, 1.6))
    a.set_ylim((0, 0.2))
chabund_ax[0].set_yticks(np.arange(0, 0.16, 0.05))
theylabs = np.round(chabund_ax[0].get_yticks() * 100)
chabund_ax[0].set_yticklabels(theylabs, fontsize = base_txtsize - 2)
chabund_ax[0].set_ylabel('Percent Channel %', fontsize = base_txtsize - 2)
chabund_ax[0].spines['left'].set_visible(True)


f.subplots_adjust(left = 0.075, right = 0.975, top = 0.925, bottom = 0.15)
# f.subplots_adjust(left = 0.05, right = 0.975)

sft = 0.02

box = chabund_ax[0].get_position()
box.x0 = box.x0 + 0.025
box.x1 = box.x1 + 0.025
chabund_ax[0].set_position(box)

box = chabund_ax[1].get_position()
box.x0 = box.x0 + 0.015
box.x1 = box.x1 + 0.015
chabund_ax[1].set_position(box)

box = xsec.get_position()
box.x0 = box.x0 + 0.01
box.x1 = box.x1 + 0.01
box.y0 = box.y0 + 0.05
box.y1 = box.y1 + 0.05
xsec.set_position(box)
xsec.set_title('(c) Example Stratigraphic Cross-section', fontsize = base_txtsize - 1, pad = 5, y = 1)

xsec.plot(0.6, 0.05, 'r.', markersize = 3)
xsec.plot(0.7, 0.045, 'r.', markersize = 3)

xsec.annotate('', xy=(0.6, 0.05),  xycoords='data',
            xytext=(0.30, 0.60), textcoords='figure fraction',
            fontsize = 8, color = 'red',
            arrowprops=dict(
                color='r', alpha = 1,
                arrowstyle="<-",
                connectionstyle="arc,angleA=0,angleB=215,armA=0,armB=75,rad=5"
#                 arrowstyle="-",  connectionstyle="arc3,rad=.2"
            ),
            horizontalalignment='right', verticalalignment='top',
            )

# "arc,angleA=0,angleB=90,armA=10,armB=10,rad=5"

xsec.annotate('', xy=(0.7, 0.045),  xycoords='data',
            xytext=(0.23, 0.26), textcoords='figure fraction',
            fontsize = 8, color = 'red',
            arrowprops=dict(
                color='r', alpha = 1,
                arrowstyle="<-",
                connectionstyle="arc,angleA=48,angleB=-90,armA=85,armB=50,rad=5"
            ),
            horizontalalignment='right', verticalalignment='top',
            )

# , bottom = 0.075, top = 0.95
f.savefig(
    fname = os.path.join(dir_path, 'analysis/chapter_2/figures/output/channel_grain_size.pdf'),
    transparent = False, dpi = 200
)

agu_data.close()
