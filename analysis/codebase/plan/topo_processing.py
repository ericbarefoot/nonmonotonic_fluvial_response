## functions for lidar data processing steps.
## Eric Barefoot
## Sep 2019

import topo_mask_gen

import os
import numpy as np
import pandas as pd
import skimage.morphology as mo
import skimage.measure as ms
from skimage.restoration import inpaint
from scipy.stats import linregress
import hdf5storage as hdf
from tqdm import tqdm

def fill_topo(topo, mask, oceanz, pit_thresh = 0.05, island_size = 5000, hole_size = 200, view_depth = 0.01, pitsize = 3):

    # we are interested in getting a mask of where we have holes in the data
    # so make a mask where everywhere that does not have data where we would want it is equal to 1

    topo_mask = np.ones(np.shape(topo))
    topo_mask[topo == 0] = 0 # mask all pixels that are zero
    topo_mask[topo <= oceanz - view_depth] = 0 # mask all areas that below a certain region below sea level
    topo_mask[~mask] = 0 # mask based on a priori decisions

    topo_mask = mo.remove_small_objects(topo_mask.astype(np.bool), min_size = island_size)
    topo_mask = mo.remove_small_holes(topo_mask, area_threshold = hole_size)
    topo_mask = mo.binary_closing(topo_mask, selem = mo.disk(10))
    topo_mask[~mask] = 0 # make the a priori mask equal to nan

    masked_topo = np.copy(topo)
    masked_topo[~topo_mask] = np.nan
    topmean = np.nanmean(masked_topo)
    masked_topo[~topo_mask] = topmean

    topo_black_hat = mo.black_tophat(masked_topo, mo.disk(pitsize))
    bhat_mask = topo_black_hat > pit_thresh
    bhat_mask[masked_topo == 0] = True

    topo_fill = inpaint.inpaint_biharmonic(masked_topo, bhat_mask)
    topo_fill[~topo_mask] = np.nan

    return topo_fill, topo_mask.astype('bool')

def get_dry_mask(topo, mask, oceanz, pit_thresh = 0.05, island_size = 5000, hole_size = 200):

    # we are interested in getting a mask of where we have holes in the data
    # so make a mask where everywhere that does not have data where we would want it is equal to 1

    dry_mask = np.ones(np.shape(topo))
    dry_mask[topo == 0] = 0 # mask all pixels that are zero
    dry_mask[topo <= oceanz] = 0 # mask all areas that below a certain region below sea level
    dry_mask[~mask] = 0 # mask based on a priori decisions

    dry_mask = mo.remove_small_objects(dry_mask.astype(np.bool), min_size = island_size)
    dry_mask = mo.remove_small_holes(dry_mask, area_threshold = hole_size)
    dry_mask = mo.binary_closing(dry_mask, selem = mo.disk(10))
    dry_mask[~mask] = 0 # make the a priori mask equal to nan

    return dry_mask.astype('bool')

def remove_mats(topo, ocean_z, mask, tol = 0.005, island_size = 4000):

    # mat_mask = np.copy(topo)
    # mat_mask = np.logical_and(topo >= ocean_z - tol, topo <= ocean_z + tol)
    mat_mask = topo >= (ocean_z - tol)
    mat_mask_no_islands = mo.remove_small_objects(mat_mask, min_size = island_size)

    mats_only = np.logical_xor(mat_mask, mat_mask_no_islands)

    mats_only = mo.binary_dilation(mats_only, selem = mo.disk(radius = 5))

    mats_only[~mo.binary_erosion(mask, selem = mo.disk(radius = 10))] = 0

    topo_mod = inpaint.inpaint_biharmonic(topo, mats_only)

    # do some additional stuff to clean up ocean bed surfaces.

    is_not_seabed = topo_mod != 0
    is_not_seabed[~mask] = 1
    not_seabed = mo.remove_small_objects(is_not_seabed, min_size = island_size)
    not_seabed = mo.binary_erosion(not_seabed, selem = mo.disk(7))
    topo_mod[~not_seabed] = 0
    topo_mod[~mask] = np.nan

    return topo_mod, mats_only

def get_slopes(calib_scans, mask, ocean_levels):

    just_below = np.zeros(np.shape(calib_scans))

    calib_scans[~mask] = 0
    calib_scans[calib_scans >= 0.65] = 0

    for i in np.arange(calib_scans.shape[0]):
        slc = fill_topo(calib_scans[i,:,:], mask)
        slc[slc > ocean_levels[i]] = np.nan
        just_below[i,:,:] = slc

    just_below[just_below == 0] = np.nan

    slopes = np.zeros(np.shape(calib_scans)[1:])

    for i in tqdm(np.arange(just_below.shape[1]), desc='calculating refraction slopes...'):
        for j in np.arange(just_below.shape[2]):
            the_slice = just_below[:,i,j]
            nonNAN = ~np.isnan(the_slice)
            if np.any(nonNAN):
                reg = linregress(ocean_levels[nonNAN] - the_slice[nonNAN], ocean_levels[nonNAN])
                slopes[i,j] = -reg[0]
            else:
                slopes[i,j] = np.nan

    slopes[np.logical_or(slopes > 3, slopes < -2)] = np.nan

    return slopes

def generate_products(refrac_calib_scans_file, ocean_z_refrac_calib_file, run_log_files, save_dir, _experiment):

    # get ocean levels hourly for whole experimental run.

    def tdb19_fun():
        _tab = pd.read_csv(run_log_file, header = 2)
        _hours = np.arange(0,987) / 24
        # get ocean levels for the postsubs refraction scans.
        # has no header in file, so supply here.
        cols = ["DateTime", "FileName", "Date", "Time", "RunTime", "Pause", "OcnZ", "TargetOcnZ", "WeirZ", "Qsed", "Qin", 'Qin On/Off', 'Qaux On/Off']
        # read in the csv
        ddat = pd.read_csv(ocean_z_refrac_calib_file, names = cols)
        # extract the right column, and scale to meters
        ocean_levels_calib = ddat['OcnZ'].values / 1000
        # specify the file path for this data numpy array.
        ocean_level_calib_file = os.path.join(save_dir, 'ocean_level_array_postsubs.npy')
        # save it
        np.save(ocean_level_calib_file, ocean_levels_calib)
        # need to get the scans from the calibration
        refrac_calib_scans = np.copy(h5py.File(refrac_calib_scans_file, 'r')['ZP'])
        dims = refrac_calib_scans.shape
        # make a mask a priori
        mask1 = topo_mask_gen.make_topo_mask1(dims, _experiment)
        # find the refraction correction values
        # _refrac_corrections = get_slopes(refrac_calib_scans, mask1, ocean_levels_calib)

        _refrac_corrections = -1 * np.ones(dims)

        _refrac_corrections[~mask1] = np.nan

        return _tab, _hours, _refrac_corrections

    def tdb15_fun():
        # raise Exception('Called ocean data subroutine for 15-2')
        names = pd.read_csv(run_log_files, header = 2, nrows = 2)
        _tab = pd.read_csv(run_log_files, skiprows = 18169)
        _tab.columns = names.columns
        _hours = (np.arange(535.4,734, 1.1) - 300) / 24

        refrac_calib_scans = hdf.loadmat(refrac_calib_scans_file)['ZD']
        dims = refrac_calib_scans.shape
        # make a mask a priori
        mask1 = topo_mask_gen.make_topo_mask1(dims, _experiment)

        _refrac_corrections = -1 * np.ones(dims[1:])

        _refrac_corrections[~mask1] = np.nan

        return _tab, _hours, _refrac_corrections

    def tdb12_fun():
        _tab = pd.read_csv(run_log_files, header = 2)
        rt = _tab['RunTime (days)']
        drt = np.diff(rt)
        pts = np.where(drt < -5)
        drt[pts[0]] = 0
        rt_fix = np.append(0,np.cumsum(drt))
        _tab['RunTime (days)'] = rt_fix

        _hours = np.arange(0,1285) / 24

        refrac_calib_scans = hdf.loadmat(refrac_calib_scans_file)['ZD']
        dims = refrac_calib_scans.shape
        # make a mask a priori
        mask1 = topo_mask_gen.make_topo_mask1(dims, _experiment)

        _refrac_corrections = -1 * np.ones(dims)

        _refrac_corrections[~mask1] = np.nan

        return _tab, _hours, _refrac_corrections
        # raise Exception('Called ocean data subroutine for 12-1')

    def getOcnData(argument):
        experiments = {
            'tdb19': tdb19_fun,
            'tdb15': tdb15_fun,
            'tdb12': tdb12_fun
        }
        ocnFunction = experiments.get(argument, 'nothing')
        return ocnFunction()

    tab, hours, refrac_corrections = getOcnData(_experiment)

    OcnZ = np.zeros(np.shape(hours))

    j = 0

    for i in hours:
        min_list = abs(tab['RunTime (days)'] - i)
        min_ind = min_list.idxmin()
        OcnZ[j] = tab['OcnZ'].iloc[min_ind]
        j += 1

    OcnZ = OcnZ / 1000

    ocean_level_file = os.path.join(save_dir, 'ocean_level_array.npy')

    np.save(ocean_level_file, OcnZ)

    # save the refraction correction arrays
    refrac_correction_file = os.path.join(save_dir, 'refraction_corrections.npy')
    np.save(refrac_correction_file, refrac_corrections)

def correct_refraction(scan_filled, ocean_level, refract_corr):

    scan_above = scan_filled >= ocean_level

    scan_filled_below = np.copy(scan_filled)
    scan_filled_below[scan_above] = np.nan
    scan_filled_below[scan_filled_below == 0] = np.nan

    scan_refract_corr = ocean_level + refract_corr * (ocean_level - scan_filled_below)

    correctSlice = np.isnan(scan_refract_corr)

    scan_corr = np.copy(scan_filled)
    scan_corr[~correctSlice] = scan_refract_corr[~correctSlice]

    return scan_corr

def swap_spare_maps(primary_array, spare_array):
    has_map = np.ones(primary_array.shape[0])
    for i in np.arange(primary_array.shape[0]):
        if np.all(primary_array[i,:,:] == 0) and not np.all(spare_array[i,:,:] == 0):
            primary_array[i,:,:] = spare_array[i,:,:]
        elif np.all(primary_array[i,:,:] == 0) and np.all(spare_array[i,:,:] == 0):
            has_map[i] = 0
    for i in np.arange(primary_array.shape[0]):
        if has_map[i] == 0:
            ii, jj = closest_non_zeros(has_map, i)
            primary_array[i,:,:] = np.mean([primary_array[ii,:,:], primary_array[jj,:,:]], axis = 0)
    return primary_array

def closest_non_zeros(l, i):
    if l[i] > 0:
        return l[i]
    arr = np.array(l)
    non_zeros = np.nonzero(arr)[0]
    distances = non_zeros - i
    neg_dist = np.float16(distances)
    neg_dist[distances > 0] = float('-inf')
    pos_dist = np.float16(distances)
    pos_dist[distances < 0] = float('inf')
    closest_idx_neg = np.min(np.where(neg_dist == np.max(neg_dist)))
    closest_idx_pos = np.min(np.where(pos_dist == np.min(pos_dist)))
    return (non_zeros[closest_idx_neg], non_zeros[closest_idx_pos])
