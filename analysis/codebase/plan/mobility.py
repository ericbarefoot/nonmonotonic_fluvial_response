# functions to calculate instantaneous mobility a la Wickert 2013
# Eric Barefoot
# Nov 2019

import numpy as np
from tqdm import tqdm

# method number 1
# scaled difference.


def overlap(K1, K2):
    w, h = K1.shape
    A = w * h
    D = np.sum(np.bitwise_xor(K1, K2))
    fw1 = np.sum(K1) / A
    fw2 = np.sum(K2) / A
    fd1 = 1 - fw1
    fd2 = 1 - fw2
    Phi = (fw1 * fd2 + fw2 * fd1)
    Ophi = 1 - (D / (A * Phi))

    return Ophi


def many_overlap(Ks, axis = 0, leavePbar = True, disablePbar = False):

    lt = Ks.shape[axis] - 1

    ov = np.full((lt, lt), np.nan)
    t0 = np.full((lt, lt), np.nan)
    t1 = np.full((lt, lt), np.nan)

    for j in tqdm(
        range(lt), desc = 'calculating overlap...',
        leave = leavePbar, disable = disablePbar
    ):
        for i in range(lt - j):
            ov[i, j] = overlap(Ks[j, ...], Ks[j + i + 1, ...])
            t0[i, j] = j
            t1[i, j] = j + i + 1

    eitherZ = np.isnan(t0) | np.isnan(t1) | np.isnan(ov)

    tt0 = t0[~eitherZ].flatten()
    tt1 = t1[~eitherZ].flatten()
    ovv = ov[~eitherZ].flatten()

    return tt0.astype(int), tt1.astype(int), ovv

# this is to be used as comparing every channel to all the channel configs that follow. then fitting a curve.

# def overlap_fit():

# Method number 2

def reworking(NKs):

    # takes masks stacked with timestep as fastest-varying index

    Fprime = np.logical_not(NKs[0])
    Aprime = np.sum(Fprime)

    nksNotReworked = NKs.cumsum(axis = 0) == 0

    Nprime = nksNotReworked.sum(axis=1).sum(axis=1)
    Freworked = 1 - (Nprime / NKs[0].size)

    return Freworked

# is a little more involved, as it involves the fluvial surface area. later...

# Method number 3


def inst_mobil(K1, K2, K3, dt = 1):

    w, h = K1.shape
    A = w * h

    D1 = np.sum(abs(K1 - K2))
    D2 = np.sum(abs(K1 - K3))

    DD = D2 - D1

    zeta = DD / (2 * A * dt)

    return zeta


def many_mobil(Ks, axis = 0):

    lt = Ks.shape[axis] - 2

    zeta = np.full(lt, np.nan)
    t0 = np.full(lt, np.nan)

    for j in tqdm(range(lt)):
        zeta[j] = inst_mobil(Ks[j, ...], Ks[j + 1, ...], Ks[j + 2, ...])
        t0[j] = j

    nanZ = np.isnan(zeta)

    return t0[~nanZ], zeta[~nanZ]
