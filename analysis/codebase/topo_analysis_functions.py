# a set of functions for AGU and preliminary paper analysis
# Eric Barefoot
# Nov 2020

import os
import numpy as np
from scipy import signal as sp
import pandas as pd
import h5py as h
# from skimage import transform as tr
from tqdm import tqdm

dir_path = os.getcwd()

# import sys
# sys.path.insert(0, dir_path)

import codebase.cube.slicedice as cut
from codebase.cube.slicedice import shallowline, steepline


# meta_data_19 = pd.read_csv(os.path.join(dir_path, 'analysis/agu_analysis/data/tdb19_metadata.csv'))
# meta_data_12 = pd.read_csv(os.path.join(dir_path, 'analysis/agu_analysis/data/tdb12_metadata.csv'))

# agu_data = h.File(os.path.join(dir_path, 'analysis/agu_analysis/data/agu_data.h5'), 'a')


def get_line(_start, _end):
    x0, y0 = _start
    x1, y1 = _end
    if abs(y1 - y0) < abs(x1 - x0):
        if x0 > x1:
            exx, why = shallowline(x1, y1, x0, y0)
        else:
            exx, why = shallowline(x0, y0, x1, y1)
    else:
        if y0 > y1:
            exx, why = steepline(x1, y1, x0, y0)
        else:
            exx, why = steepline(x0, y0, x1, y1)
    dx = np.diff(exx)
    dy = np.diff(why)
    d = np.insert(np.sqrt(dx * dx + dy * dy), 0, 0).cumsum()
    return exx, why, d


def fin(theta, ymax, xmax, cc):
    radtheta = np.radians(theta)
    slp = np.tan(radtheta)
    if (theta > 0) & (theta < 90):
        rmax = np.sqrt((xmax - cc[0]) ** 2 + (ymax - cc[1]) ** 2)
        xproj = np.cos(radtheta) * rmax + cc[0]
        yproj = np.sin(radtheta) * rmax + cc[1]
        if (yproj < ymax) & (xproj > xmax):
            y = int(yproj)
            x = int(xmax)
        elif (yproj > ymax) & (xproj < xmax):
            y = int(ymax)
            x = int(xproj)
    elif (theta < 0) & (theta > -90):
        ytst = xmax * slp + cc[1]
        if (ytst > 0):
            y = int(ytst)
            x = int(xmax)
        elif (ytst < 0):
            y = 0
            x = int(-(1/slp) * cc[1] + cc[0])
    elif (theta > 90) & (theta < 180):
        xtst = ymax * 1/slp + cc[0]
        if (xtst > 0):
            y = int(ymax)
            x = int(xtst)
        elif (xtst < 0):
            y = int(-(slp * cc[0]) + cc[1])
            x = 0
    elif theta == 0:
        x = xmax
        y = cc[1]
    elif theta == 90:
        x = cc[0]
        y = ymax
    if x > xmax:
        x = xmax
    if y > ymax:
        y = ymax
    return x, y


def detrend_cone(surf, _apex, corner1, corner2, dtheta, _disable_pbar = False):
    ''' subtract a deformable conical surface from the delta '''
    tdest = np.copy(surf)
    cx, cy = _apex
    cx1, cy1 = corner1
    cx2, cy2 = corner2
    thi = np.degrees(np.arctan((cy1 - cy) / (cx1 - cx)))
    thf = np.degrees(np.arctan((cy2 - cy) / (cx2 - cx))) + 180
    angles = np.arange(thi-1, thf+1, dtheta)
    for i in tqdm(angles, leave = False, disable = _disable_pbar):
        ex, ey = fin(i, surf.shape[0], surf.shape[1], (cx,cy))
        x, y, d = get_line((cx,cy), (ex,ey))
        sig = surf[y, x]
        nn = ~np.isnan(sig)
        if np.any(nn):
            sigg = sig[nn]
            dt = sp.detrend(sigg)
            tdest[y[nn], x[nn]] = dt
    return tdest


def detrend_average_cone(surf, _apex, corner1, corner2, ntheta, dtheta, nanthresh = 0.2):
    ''' subtract an average conical surface from the delta '''
    tdest = np.copy(surf)
    cx, cy = _apex
    cx1, cy1 = corner1
    cx2, cy2 = corner2
    thi = np.degrees(np.arctan((cy1 - cy) / (cx1 - cx)))
    thf = np.degrees(np.arctan((cy2 - cy) / (cx2 - cx))) + 180
    angles = np.linspace(thi-1, thf+1, ntheta)
    slopes = np.zeros(angles.shape)
    intercepts = np.zeros(angles.shape)
    for j, i in enumerate(tqdm(angles, disable=(len(angles)<5000))):
        ex, ey = fin(i, surf.shape[0], surf.shape[1], (cx,cy))
        x, y, d = get_line((cx,cy), (ex,ey))
        sig = surf[y, x]
        nn = ~np.isnan(sig)
        if np.sum(nn) > nanthresh * len(d):
            slopes[j], intercepts[j] = np.polyfit(d[nn], sig[nn], 1)
    mean_slope = slopes.mean()
    mean_intercept = intercepts.mean()
    angles = np.arange(thi-1, thf+1, dtheta)
    for i in tqdm(angles):
        ex, ey = fin(i, surf.shape[0], surf.shape[1], (cx,cy))
        x, y, d = get_line((cx,cy), (ex,ey))
        sig = surf[y, x]
        nn = ~np.isnan(sig)
        if np.sum(nn) > nanthresh * len(d):
            tdest[y[nn], x[nn]] = sig[nn] - (mean_intercept + d[nn] * mean_slope)
    # for j, i in enumerate(tqdm(angles, disable=(len(angles)<5000))):
    #     ex, ey = fin(i, surf.shape[0], surf.shape[1], (cx,cy))
    #     x, y, d = get_line((cx,cy), (ex,ey))
    #     sig = surf[y, x]
    #     nn = ~np.isnan(sig)
    #     if np.sum(nn) > nanthresh * len(d):
    #         tdest[y[nn], x[nn]] = mean_intercept + d[nn] * mean_slope
    return tdest


# measure conical slope

def measure_slopes(surf, _apex, corner1, corner2, ntheta, nanthresh = 0.2):
    ''' measure the slope of a conical surface from the delta '''
    cx, cy = _apex
    cx1, cy1 = corner1
    cx2, cy2 = corner2
    thi = np.degrees(np.arctan((cy1 - cy) / (cx1 - cx)))
    thf = np.degrees(np.arctan((cy2 - cy) / (cx2 - cx))) + 180
    angles = np.linspace(thi-1, thf+1, ntheta)
    slopes = np.zeros(angles.shape)
    for j, i in enumerate(tqdm(angles, disable=(len(angles)<5000))):
        ex, ey = fin(i, surf.shape[0], surf.shape[1], (cx,cy))
        x, y, d = get_line((cx,cy), (ex,ey))
        sig = surf[y, x]
        nn = ~np.isnan(sig)
        if np.sum(nn) > nanthresh * len(d):
            slopes[j] = np.polyfit(d[nn]/200, sig[nn], 1)[0]
    return slopes, angles


def get_several_radii(surf, _apex, corner1, corner2, ntheta, nanthresh = 0.2):
    ''' return line paths for several '''
    cx, cy = _apex
    cx1, cy1 = corner1
    cx2, cy2 = corner2
    thi = np.degrees(np.arctan((cy1 - cy) / (cx1 - cx)))
    thf = np.degrees(np.arctan((cy2 - cy) / (cx2 - cx))) + 180
    angles = np.linspace(thi-1, thf+1, ntheta)
    slopes = np.zeros(angles.shape)
    for j, i in enumerate(tqdm(angles, disable=(len(angles)<5000))):
        ex, ey = fin(i, surf.shape[0], surf.shape[1], (cx,cy))
        x, y, d = get_line((cx,cy), (ex,ey))
        sig = surf[y, x]
        nn = ~np.isnan(sig)
        if np.sum(nn) > nanthresh * len(d):
            slopes[j] = np.polyfit(d[nn], sig[nn], 1)[0]
    return slopes, angles

# measure relief

# def

# measure shoreline rugosity

# def
