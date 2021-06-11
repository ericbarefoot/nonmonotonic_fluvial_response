##
## functions to calculate compensation statistics in 1, 2, and 3D
## EA Barefoot
## Nov 2019

import numpy as np
from tqdm import tqdm

def comp_3D(array, st=1, nobar = False):
    # takes a 3D numpy array where values are z values, first axis is time,
    # and second two axes are spatial coordinates.
    n = array.shape[0]
    if np.any(np.isnan(array)):
        def mymean(x):
            return np.nanmean(x)
        def mystd(x):
            return np.nanstd(x)
    else:
        def mymean(x):
            return np.mean(x)
        def mystd(x):
            return np.std(x)
    T = mymean(array[-1,:] - array[0,:])
    siglength = np.int((n-1) * (1 + (n-1)) / 2)
    sig = np.zeros(siglength)
    eta = np.zeros(siglength)
    k = 0
    pbar = tqdm(total = siglength, disable = nobar)
    for i in np.arange(1,n,st):
        for j in np.arange(0,i,st):
            deta = array[i,:,:] - array[j,:,:]
            detabar = T/n * (i-j)
            zeta = deta / detabar
            eta[k] = detabar
            sig[k] = mystd(zeta)
            k += 1
            pbar.update(1)
    pbar.close()
    return sig, eta

def comp_2D(array, st = 1, nobar = False):
    # takes a 2D numpy array where values are z values
    # first axis are t coordinates, and second axis are x coordinates
    n = array.shape[0]
    if np.any(np.isnan(array)):
        def mymean(x):
            return np.nanmean(x)
        def mystd(x):
            return np.nanstd(x)
    else:
        def mymean(x):
            return np.mean(x)
        def mystd(x):
            return np.std(x)
    T = mymean(array[-1,:] - array[0,:])
    siglength = np.int((n-1) * (1 + (n-1)) / 2)
    sig = np.zeros(siglength)
    eta = np.zeros(siglength)
    k = 0
    pbar = tqdm(total = siglength, disable = nobar)
    for i in np.arange(1,n):
        for j in np.arange(0,i):
            deta = array[i,:] - array[j,:]
            detabar = T/n * (i-j)
            zeta = deta / detabar
            # zzeta = zeta * zeta
            # intZeta = np.trapz(y = zzeta, x = dz)
            eta[k] = detabar
            # sig[k] = np.sqrt(intZeta)
            sig[k] = mystd(zeta)
            k += 1
            pbar.update(1)
    pbar.close()
    return sig, eta

def comp_1D(array, st = 1):
    # takes a 1D numpy array where values are z values
    #
    # T = np.max(array, axis = 0) - np.min(array, axis = 0)
    n = array.shape[-1]
    #
    siglength = np.int((n-1) * (1 + (n-1)) / 2)
    deta = np.zeros(siglength)
    zeta = np.zeros(siglength)
    dif = np.zeros(siglength)
    k = 0
    pbar = tqdm(total = siglength + n)
    for i in np.arange(1,n):
        for j in np.arange(0,i):
            # detabar = (T / n) * (i - j)
            deta[k] = array[i] - array[j]
            # zeta[k] = deta[k] / detabar
            dif[k] = i - j
            k += 1
            pbar.update(1)
    dat = np.vstack((dif, deta)).T
    dif_u = np.unique(dif, return_counts=True)
    enough = dif_u[1] > 5
    dif_unique = dif_u[0]
    dif_l = len(dif_unique[enough])
    sig = np.zeros(dif_l)
    eta = np.zeros(dif_l)
    k = 0
    for i in dif_unique[enough]:
        sl = dat[dat[:,0] == i]
        eta[k] = np.mean(sl[:,1])
        sig[k] = np.std(sl[:,1]) / eta[k]
        k += 1
        pbar.update(1)
    pbar.close()
    return sig, eta
