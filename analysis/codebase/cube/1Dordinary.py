
import numpy as np
from scipy import signal as sig
import matplotlib.pyplot as plt
import matplotlib
import codebase.cube.slicedice as cut
matplotlib.use('Qt4Agg')
plt.ion()

eta = np.load('/Users/ericfoot/Dropbox/research/transition_time_experiments/data/derived_data/tdb19/python_arrays/processed_scans.npy')
syn = np.load('/Users/ericfoot/Dropbox/research/transition_time_experiments/data/derived_data/tdb19/python_arrays/synstrat_tdb19_20200414.npy')

eta12 = np.load('/Users/ericfoot/Dropbox/research/transition_time_experiments/data/derived_data/tdb12/python_arrays/processed_scans.npy')
syn12 = np.load('/Users/ericfoot/Dropbox/research/transition_time_experiments/data/derived_data/tdb12/python_arrays/synstrat_tdb12_20200415.npy')

def make_strat_column():
    _eta = np.random.normal(size = 10000).cumsum()
    _eta = _eta - np.mean(_eta)
    _syn = np.flip(np.minimum.accumulate(np.flip(_eta)))
    _bool = _eta == _syn
    return _eta, _syn, _bool

def get_og_percentiles(_eta):
    sp.signal.detrend()
    _eta = _eta - np.mean(_eta)
    _syn = np.flip(np.minimum.accumulate(np.flip(_eta)))
    _bool = _eta == _syn
    return (_eta.argsort().argsort() + 1) / _eta.size

def ordinariness1D(_per, _bool):
    preserved = _per[_bool]
    return 1 - 2 * np.median(preserved)

coreX = np.random.randint(100, 200, 5)
coreY = np.random.randint(100, 200, 5)

plt.cla()
plt.imshow(syn[300])
plt.plot(coreX, coreY, 'ro')

plt.cla()
plt.plot(eta[...,coreX,coreY])
plt.plot(sig.detrend(eta[...,coreX[0],coreY[0]]))
