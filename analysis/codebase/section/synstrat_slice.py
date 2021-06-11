# a simple script to make a nice plot of an arbitrary slice

# import analysis.TDB_topo_processing.python_scripts.
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
# matplotlib.use('Qt4Agg')
import os

def sliceplot(_synstrat, _slice, xx, yy, path, file, map_i = 1, highlight_int = False, show_map_line = False, indVec = None):
    nan_cols = np.all(np.isnan(_slice), axis = 0)
    ss = _slice[:,~nan_cols]

    fig3 = plt.figure(constrained_layout=True)
    gs = fig3.add_gridspec(1, 3)
    f3_ax1 = fig3.add_subplot(gs[0,0])
    f3_ax1.set_title('Map View')
    f3_ax2 = fig3.add_subplot(gs[0,1:])
    f3_ax2.set_title('Synthetic Stratigraphy')

    mapimage = f3_ax1.imshow(_synstrat[map_i,...], cmap = 'Greys_r')
    refline = f3_ax1.plot(yy[~nan_cols], xx[~nan_cols], 'r')
    slcplot = f3_ax2.plot(np.transpose(ss), 'k', alpha = 1, linewidth = 0.3)
    if highlight_int:
        highlightplot = f3_ax2.plot(np.transpose(ss[indVec,...]), 'r', alpha = 0.4, linewidth = 0.3)

    if show_map_line:
        mapIplot = f3_ax2.plot(np.transpose(ss[map_i,...]), 'b', alpha = 1, linewidth = 1)

    fig3.set_size_inches(15, 5)

    plt.savefig(os.path.join(path, file))
