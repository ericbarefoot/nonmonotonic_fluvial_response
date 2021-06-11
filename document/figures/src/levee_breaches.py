import os
import numpy as np
import pandas as pd
import h5py as h
from tqdm import tqdm
import seaborn as sns

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
# agu_data.visit(print)


i = 26
tind = agu_data['topoqv15'].attrs['IDs'][i]
iind = int(np.where(agu_data['imagRefqv15'].attrs['IDs'] == (tind - 4))[0])
mind = int(np.where(agu_data['channelMasksqv15'].attrs['IDs'] == (tind - 4))[0])
mmind = int(np.where(agu_data['manualChannelMasksqv15'].attrs['IDs'] == (tind - 4))[0])
img = agu_data['imagRefqv15'][iind]
# tind
allbrchData = pd.read_csv(os.path.join(dir_path, 'data/raw_data/breach_counting_all.csv'))
thisbreach = allbrchData[allbrchData['img'] == tind + 3]
# print(tind - 2)
# print(tind+3)
# print(allbrchData['img'].unique())


dat = pd.read_csv(os.path.join(dir_path, 'data/derived_data/breach_counting_counts.csv'))
sns.set_style("ticks")

f = plt.figure(figsize = (6,3.5))

gg = f.add_gridspec(ncols = 2, nrows = 2,
                    height_ratios = [1,1], width_ratios = [1,3],
#                     wspace=0.0, hspace=0.1
                   )
base_fontsize = 10

image = f.add_subplot(gg[:, 1])
image.spines['top'].set_visible(False)
image.spines['bottom'].set_visible(False)
image.spines['right'].set_visible(False)
image.spines['left'].set_visible(False)
# image.set_title('Image', fontsize = base_fontsize)
image.set_xticks(np.arange(0, 2.6, 0.5))
image.set_xticklabels(image.get_xticks(), fontsize = base_fontsize - 2)
image.set_xlabel('Distance (m)', fontsize = base_fontsize - 1)
image.set_yticks(np.arange(0, 2.6, 0.5))
image.set_yticklabels(image.get_yticks(), fontsize = base_fontsize - 2)
image.set_ylabel('Distance (m)', fontsize = base_fontsize - 1)
image.imshow(img,
            extent = [0, (img.shape[1] / 200), 0, (img.shape[0] / 200)],
            origin = 'lower'
          )

image.plot(thisbreach['y'], thisbreach['x'], '.', color = '#ff0000', markersize = 3)

brchs = f.add_subplot(gg[0,0])

brchs = sns.pointplot(
    x="qv", y="n", data=dat,
    color = 'k', ax = brchs,
    scale = 0.75, errwidth = 2
)
sns.despine()

brchs.set_ylabel('Number of\nLevee Breaches', fontsize = base_fontsize- 1)
brchs.set_xlabel('Flood Intensity $(Q_v)$', fontsize = base_fontsize- 1)
# brchs.set_xticks(np.arange(1, 3.1, 1))
# brchs.set_xticklabels(brchs.get_xticks(), size = base_fontsize - 2)
# brchs.set_xlabel('Distance (m)', fontsize = base_fontsize - 1)
brchs.set_yticks(np.arange(4, 12.1, 1))
brchs.set_yticklabels(brchs.get_yticks(), fontsize = base_fontsize - 2)
# brchs.set_ylabel('Distance (m)', fontsize = base_fontsize - 1)
sns.set(font_scale = 0.75)

image.annotate('levee breaches\nwith water escaping\nat low flow',
            xy=(0.734, 0.898),  xycoords='data',
            xytext=(0.30, 0.35), textcoords='figure fraction',
            fontsize = 8, color = '#ff0000',
            arrowprops=dict(
                color='#ff0000', alpha = 1,
                arrowstyle="-",  connectionstyle="arc3,rad=0"
                # arrowstyle="<-",  connectionstyle="arc,angleA=45,angleB=215,armA=97,armB=118,rad=5"
            ),
            horizontalalignment='right', verticalalignment='top',
            )

image.annotate('levees with annealed\nbreaches',
            xy=(0.55, 0.45),  xycoords='data',
            xytext=(0.30, 0.15), textcoords='figure fraction',
            fontsize = 8, color = '#ff0000',
            arrowprops=dict(
                color='#ff0000', alpha = 1, relpos = (1, 0.5),
                arrowstyle="-",  connectionstyle="arc3,rad=0"
                # arrowstyle="<-",  connectionstyle="arc,angleA=45,angleB=215,armA=97,armB=118,rad=5"
            ),
            horizontalalignment='right', verticalalignment='top',
            )

f.subplots_adjust(top = 0.925, right = 1.05)

f.savefig(fname = os.path.join(dir_path, 'analysis/chapter_2/figures/output/breaches.pdf'),
          transparent = False, dpi = 200)
