# script to extract compensation measurements for sections of the synthetic
# stratigraphy and export to array, then read in to R and calculate the
# Eric Barefoot
# Feb 2020

import codebase.cube.compensation as cp
import codebase.cube.slicedice as cut

import numpy as np
import pandas as pd
import scipy.special as sp
import os
import argparse
from tqdm import tqdm
# want there to be option to do radial slices of increasing sizes and also 'cores'
# of arbitrary size and location.
# for now, just start with radial slices of increasing size.

# TODO: make the function as it is now about extracting at radii, then
# make another one that extracts for a given time interval. This one can be specific to tdb19 for now.

def group(lst, n):
    for i in range(0, len(lst), n):
        val = lst[i:i+n]
        if len(val) == n:
            yield tuple(val)

def calcCompAtRadii(_synstrat, **kwargs):

    _apex = kwargs.pop('apex')
    Nsurfaces = _synstrat.shape[0]
    Npairs_per = int(sp.comb(Nsurfaces, 2))

    if 'radii' in kwargs:
        radii = kwargs.pop('radii')
    elif 'radii_seq' in kwargs:
        _radii = kwargs.pop('radii_seq')
        radii = np.arange(_radii[0], _radii[1], _radii[2])

    Npairs = len(radii) * Npairs_per

    sig = np.zeros(Npairs)
    eta = np.zeros(Npairs)
    rad = np.zeros(Npairs)
    j = 0
    for i in radii:
        slc, xx, yy = cut.circular_slice(_synstrat, _apex, i)
        nan_cols = np.all(np.isnan(slc), axis = 0)
        ss = slc[:,~nan_cols]
        sigTemp, etaTemp = cp.comp_2D(ss[:,2:-2])
        sig[j:(j+Npairs_per)] = sigTemp
        eta[j:(j+Npairs_per)] = etaTemp
        rad[j:(j+Npairs_per)] = i
        j += Npairs_per

    dat = pd.DataFrame(np.vstack((sig,eta,rad)).T)
    dat.columns = ['sigma', 'eta', 'radius']

    return dat

def calcCompAtChords(_synstrat, **kwargs):
    _starts = kwargs.pop('starts')
    _ends = kwargs.pop('ends')
    starts = list(group(_starts, 2))
    ends = list(group(_ends, 2))
    Nsurfaces = _synstrat.shape[0]

    # print(zip(starts, ends))
    #
    # for x, y in zip(starts, ends):
    #     print(x)
    # return x

    Npairs_per = int(sp.comb(Nsurfaces, 2))

    Npairs = len(starts) * Npairs_per

    sig = np.zeros(Npairs)
    eta = np.zeros(Npairs)
    srt = np.zeros((Npairs,2))
    end = np.zeros((Npairs,2))
    sec = np.zeros(Npairs)
    j = 0
    k = 1
    for the_srt, the_end in zip(starts, ends):
        slc, xx, yy = cut.chord_slice(_synstrat, the_srt, the_end)
        nan_cols = np.all(np.isnan(slc), axis = 0)
        ss = slc[:,~nan_cols]
        sigTemp, etaTemp = cp.comp_2D(ss[:,2:-2])
        sig[j:(j+Npairs_per)] = sigTemp
        eta[j:(j+Npairs_per)] = etaTemp
        srt[j:(j+Npairs_per),:] = the_srt
        end[j:(j+Npairs_per),:] = the_end
        sec[j:(j+Npairs_per)] = k
        j += Npairs_per
        k += 1
    sig = np.reshape(sig,(Npairs,1))
    eta = np.reshape(eta,(Npairs,1))
    sec = np.reshape(sec,(Npairs,1))

    dat = pd.DataFrame(np.hstack((sig,eta,srt,end,sec)))
    dat.columns = ['sigma', 'eta', 'start_x', 'start_y', 'end_x', 'end_y', 'section']

    return dat

def calcCompCircForTimeInt(_synstrat, **kwargs):

    _apex = kwargs.pop('apex')
    _radius = kwargs.pop('radius')
    _times = kwargs.pop('times')
    _labels = kwargs.pop('labels')
    nobar = kwargs.pop('disable_bar')

    time_slices = list(group(_times,2))

    Nsurfaces = np.zeros(len(time_slices), dtype = 'int')
    Npairs_per = np.zeros(len(time_slices), dtype = 'int')

    for i in range(len(time_slices)):
        Nsurfaces[i] = _synstrat[np.arange(time_slices[i][0], time_slices[i][1]),:,:].shape[0]
        Npairs_per[i] = int(sp.comb(Nsurfaces[i], 2))

    # Npairs_per = Npairs_per.astype('int')
    Npairs = int(np.sum(Npairs_per))

    sig = np.zeros(Npairs)
    eta = np.zeros(Npairs)
    t_srt = np.zeros(Npairs)
    t_end = np.zeros(Npairs)
    rad = np.zeros(Npairs)
    lab = pd.DataFrame('', index = range(Npairs), columns = ['interval'])
    j = 0
    for i in range(len(time_slices)):
        slc, xx, yy = cut.circular_slice(_synstrat, _apex, _radius)
        nan_cols = np.all(np.isnan(slc), axis = 0)
        ss = slc[~nan_cols]
        st = ss[np.arange(time_slices[i][0], time_slices[i][1]),:]
        sigTemp, etaTemp = cp.comp_2D(st[:,2:-2], nobar = nobar)
        sig[j:(j+Npairs_per[i])] = sigTemp
        eta[j:(j+Npairs_per[i])] = etaTemp
        t_srt[j:(j+Npairs_per[i])] = time_slices[i][0]
        t_end[j:(j+Npairs_per[i])] = time_slices[i][1]
        rad[j:(j+Npairs_per[i])] = _radius
        lab.iloc[j:(j+Npairs_per[i]),:] = _labels[i]
        j += Npairs_per[i]

    dat = pd.DataFrame(np.vstack((sig,eta,t_srt,t_end,rad)).T)
    dat.columns = ['sigma', 'eta', 't_start', 't_end', 'radius']
    dat = dat.join(lab)

    return dat

def calcCompAtRadiiForTimeInt(_synstrat, **kwargs):

    _apex = kwargs.pop('apex')
    if 'radii' in kwargs:
        radii = kwargs.pop('radii')
    elif 'radii_seq' in kwargs:
        _radii = kwargs.pop('radii_seq')
        radii = np.arange(_radii[0], _radii[1], _radii[2])
    _times = kwargs.pop('times')
    _labels = kwargs.pop('labels')

    data = calcCompCircForTimeInt(_synstrat, apex = _apex, radius = radii[0], times = _times, labels = _labels, disable_bar = True)

    for i in tqdm(radii[1:]):
        dat = calcCompCircForTimeInt(_synstrat, apex = _apex, radius = i, times = _times, labels = _labels, disable_bar = True)
        data = data.append(dat, ignore_index = True)

    return data

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Options for computing and exporting compensational stats.')

    parser.add_argument('--root_dir', '-r', help = 'the path to the root directory of the project.')
    parser.add_argument('--experiment', '-e', choices = ['tdb12', 'tdb15', 'tdb19'], help = 'a tag for the experiment id that is being processed, one of ("tdb19", "tdb15", or "tdb12")')
    parser.add_argument('--data_file', '-d', help = 'the name of the synthetic stratigraphy data file. Makes strong assumptions about the structure of the datasets.')
    parser.add_argument('--out_file', '-o', help = 'the name of the data file to export to R. Makes strong assumptions about the structure of the datasets.')
    parser.add_argument('--function', '-f', help = 'the name of the function in here that you want to calculate.')

    parser.add_argument('--circle_center', type = int, nargs = 2)
    parser.add_argument('--radius_from_to_by', type = int, nargs = 3)
    parser.add_argument('--radii', type = int, nargs = '*')

    parser.add_argument('--start_points', type = int, nargs = '*')
    parser.add_argument('--end_points', type = int, nargs = '*')

    parser.add_argument('--times', type = int, nargs = '*')
    parser.add_argument('--time_labels', nargs = '*')
    parser.add_argument('--radius', type = int)

    args = parser.parse_args()

    func_list = [calcCompAtRadii, calcCompAtChords, calcCompCircForTimeInt, calcCompAtRadiiForTimeInt]
    func_names = [f.__name__ for f in func_list]
    funcs_dict = dict(zip(func_names, func_list))

    f = funcs_dict[args.function]

    # print(args.time_labels)
    # print(args.times)

    data_file = os.path.join(args.root_dir, 'data', 'derived_data', args.experiment, 'python_arrays', args.data_file)

    synstrat = np.load(data_file)

    cData = f(synstrat, apex = args.circle_center, radii_seq = args.radius_from_to_by, radius = args.radius, labels = args.time_labels, starts = args.start_points, ends = args.end_points, times = args.times, radii = args.radii)

    cData.to_csv(os.path.join(args.root_dir, 'data', 'derived_data', args.experiment, 'r_data', args.out_file), index = False)

    # np.save(os.path.join(args.root_dir, 'data', 'derived_data', args.experiment, 'r_data', args.out_file), cData)
