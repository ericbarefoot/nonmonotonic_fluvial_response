import os
import numpy as np
import h5py as h
import matplotlib
from matplotlib import pyplot as plt

dir_path = os.getcwd()

agu_data = h.File(os.path.join(dir_path, 'data/raw_data/agu_data.h5'), 'r')

i = 85
tind = agu_data['topoqv15'].attrs['IDs'][i]
iind = int(np.where(agu_data['imagRefqv15'].attrs['IDs'] == (tind - 4))[0])
mind = int(np.where(agu_data['channelMasksqv15'].attrs['IDs'] == (tind - 4))[0])
mmind = int(np.where(agu_data['manualChannelMasksqv15'].attrs['IDs'] == (tind - 4))[0])

f = plt.figure(figsize = (6, 3))

gg = f.add_gridspec(
    ncols = 3, nrows = 1,
    wspace=0.1, hspace=0.1
) #  width_ratios = [1,1,1]

dd = [
    agu_data['imagRefqv15'][iind],
    agu_data['channelMasksqv15'][mind],
    agu_data['manualChannelMasksqv15'][mmind]

]
tt = ['Image', 'Automatic Mask', 'Manual Retouch']

# img.set_ylabel('Distance (m)')

base_fontsize = 8

img = f.add_subplot(gg[0, 0])
aut = f.add_subplot(gg[0, 1])
msk = f.add_subplot(gg[0, 2])

ax = [img, aut, msk]

j = 0
for dif, tit, i in zip(dd, tt, ax):
    # i.set_facecolor('#000000')
    i.set_title(tit, fontsize = base_fontsize)
    i.spines['top'].set_visible(False)
    i.spines['bottom'].set_visible(False)
    i.spines['right'].set_visible(False)
    i.spines['left'].set_visible(False)
    # i.tick_params(colors='#FFFFFF')
    # i.title.set_color('#FFFFFF')
    # i.xaxis.label.set_color('#FFFFFF')
    # i.yaxis.label.set_color('#FFFFFF')
    # i.set_xlabel('Distance (m)')
    i.set_xticks(np.arange(0, 2.6, 0.5))
    i.set_xticklabels(i.get_xticks(), fontsize = base_fontsize - 2)
    i.set_xlabel('Distance (m)', fontsize = base_fontsize - 1)
    i.set_yticks([])
    if j == 0:
        i.imshow(dif,
            extent = [0, (dif.shape[1] / 200), 0, (dif.shape[0] / 200)],
            origin = 'lower')
    else:
        i.imshow(dif,
            extent = [0, (dif.shape[1] / 200), 0, (dif.shape[0] / 200)],
            origin = 'lower', cmap = 'Greys_r')
    j += 1

ax[0].set_yticks(np.arange(0, 2.6, 0.5))
ax[0].set_yticklabels(ax[0].get_yticks(), fontsize = base_fontsize - 2)
ax[0].set_ylabel('Distance (m)', fontsize = base_fontsize - 1)

f.subplots_adjust(left = 0.075, bottom = 0.075, top = 0.925, right = 0.975)

f.savefig(fname = os.path.join(dir_path, 'analysis/chapter_2/figures/output/channel_mask_procedure.pdf'), transparent = False, dpi = 200)

f.savefig(fname = os.path.join(dir_path, 'analysis/chapter_2/figures/output/channel_mask_procedure.png'), transparent = True, dpi = 200)

agu_data.close()
