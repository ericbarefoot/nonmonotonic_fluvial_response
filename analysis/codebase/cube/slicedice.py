#
# functions to separate synthetic stratigraphy into chunks and slices for analysis
# Eric Barefoot
# Feb 2020

import numpy as np
np.seterr(divide='ignore', invalid='ignore')

def circular_slice(_array, _center, _radius):
    exx, why = circle(_center, _radius)
    ind = np.logical_and(exx >= 0, why >= 0)

    arr = np.copy(_array)
    slc = arr[:, exx[ind], why[ind]]

    return slc, exx[ind], why[ind]

def chord_slice(_array, _start, _end):
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

    arr = np.copy(_array)
    slc = arr[:, exx, why]

    return slc, exx, why

def shallowline(x0, y0, x1, y1):
    dx = x1 - x0
    dy = y1 - y0
    yi = 1
    if dy < 0:
    	yi = -1
    	dy = -dy
    D = 2 * dy - dx
    x = np.arange(x0, x1, 1, dtype = 'float')
    y = np.zeros(len(x))
    yy = y0
    for i in np.arange(len(x)):
    	y[i] = yy
    	if D > 0:
    		yy = yy + yi
    		D = D - 2 * dx
    	D = D + 2 * dy
    xI = np.argsort(x)
    x = x[xI]
    y = y[xI]
    return x.astype('int'), y.astype('int')

def steepline(x0, y0, x1, y1):
    dx = x1 - x0
    dy = y1 - y0
    xi = 1
    if dx < 0:
    	xi = -1
    	dx = -dx
    D = 2 * dx - dy
    y = np.arange(y0, y1, 1, dtype = 'float')
    x = np.zeros(len(y))
    xx = x0
    for i in np.arange(len(y)):
    	x[i] = xx
    	if D > 0:
    		xx = xx + xi
    		D = D - 2 * dy
    	D = D + 2 * dx
    yI = np.argsort(y)
    x = x[yI]
    y = y[yI]
    return x.astype('int'), y.astype('int')

def circle(center, radius):

    circ = 2 * np.pi * radius

    x0, y0 = (0,0)
    f = 1 - radius
    ddf_x = 1
    ddf_y = -2 * radius
    x = 0
    y = radius

    xc = np.zeros(int(circ))
    yc = np.zeros(int(circ))

    # xc[0] = x0
    # yc[0] = y0 + radius
    xc[0] = x0
    yc[0] = y0 - radius
    xc[1] = x0 + radius
    yc[1] = y0

    j = 2

    while x < y:
        if f >= 0:
            y -= 1
            ddf_y += 2
            f += ddf_y
        x += 1
        ddf_x += 2
        f += ddf_x
        xc[j] = x0 + x
        yc[j] = y0 + y
        j += 1
        xc[j] = x0 + x
        yc[j] = y0 - y
        j += 1
        xc[j] = x0 + y
        yc[j] = y0 + x
        j += 1
        xc[j] = x0 + y
        yc[j] = y0 - x
        j += 1

    not_filled = np.logical_and(xc == 0, yc == 0)

    xc = xc[~not_filled]
    yc = yc[~not_filled]

    ycI = np.argsort(np.arctan(yc/xc))

    xc_2 = xc[ycI]
    yc_2 = yc[ycI]

    xc_1 = -xc_2[1:-1]
    yc_1 = yc_2[1:-1]

    xc = np.append(xc_2, np.flip(xc_1))
    yc = np.append(yc_2, np.flip(yc_1))

    xc = xc + center[0]
    yc = yc + center[1]

    return xc.astype('int'), yc.astype('int')
