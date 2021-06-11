## script for loading and processing lidar data from Tulane Experiments
## Eric Barefoot
## Aug 2019

# Import local libraries

import topo_processing
import topo_mask_gen

# Import system libraries

import os
import numpy as np
from tqdm import tqdm
import h5py
import argparse
from datetime import datetime

# Build file paths (TODO:supply root_dir as argument to function or command line.)
# structure of project should stay the same, so coding these relative paths is fine.

# MAJOR TODO is to make this process dataset-agnostic.
# in particular - what about specifying a command line argument for scans to exclude, etc.

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def process_lidar(root_dir, _experiment, _regenerate = False):

    raw_data_dir = os.path.join(root_dir, 'data', 'raw_data', _experiment)
    python_data_dir = os.path.join(root_dir, 'data', 'derived_data', _experiment, 'python_arrays')
    matlab_data_dir = os.path.join(root_dir, 'data', 'derived_data', _experiment, 'matlab_3d_arrays')
    fig_dir = os.path.join(root_dir, 'figures', 'quick_figs')
    wScanFile = os.path.join(matlab_data_dir, 'wet', 'ZdataW.mat')
    dScanFile = os.path.join(matlab_data_dir, 'dry', 'ZdataD.mat')
    ocean_level_file = os.path.join(python_data_dir, 'ocean_level_array.npy')
    refraction_correction_file = os.path.join(python_data_dir, 'refraction_corrections.npy')

    print('reading in scan data...')

    if os.path.exists(dScanFile):
        dryScan = np.copy(h5py.File(dScanFile)['ZD'])
    if os.path.exists(wScanFile):
        wetScan = np.copy(h5py.File(wScanFile)['ZW'])

    # logic to determine whether to re-run all code
        # i.e. whether to re-run regressions on postsubs, etc.
        # if '--generate' supplied from command line call, then re-run prep scripts. otherwise load from disk

    # This function is for TDB_19_1
    def tdb19_gen():
        print('generating necessary products...')
        _spurious_scans_dry = [274]
        _spurious_scans_wet = None
        if _regenerate:
            # list files necessary to generate data products
            ocean_z_refract_calib_file = os.path.join(raw_data_dir, 'refraction_FARO_times.csv')
            postScanFile = os.path.join(matlab_data_dir, 'ZdataP.mat')
            run_log_file = os.path.join(raw_data_dir, 'TDB_19_1_Subside_log.csv')
            # generates:
            # ocean levels for every hour
            # ocean levels during calibration scans
            # corrections for refraction.
            topo_processing.generate_products(postScanFile, ocean_z_refract_calib_file, run_log_file, python_data_dir, _experiment)
            # then load them from disk because the function returns nothing
            _ocean_levels = np.load(ocean_level_file)
            _ref_corrections = np.load(refraction_correction_file)
        else:
            # otherwise, load them from disk
            _ocean_levels = np.load(ocean_level_file)
            _ref_corrections = np.load(refraction_correction_file)
        return _ocean_levels, _ref_corrections, _spurious_scans_dry, _spurious_scans_wet
    # This function is for TDB_15_1
    def tdb15_gen():
        print('generating necessary products...')
        _spurious_scans_dry = [41,42,43,44,45,55,56]
        _spurious_scans_wet = [41,42,43,44,45,46,55,56]
        # print('called TDB_15_2 regeneration function!')
        if _regenerate:
            # raise Exception('called TDB_15_2 regeneration function!')
            # list files necessary to generate data products
            # ocean_z_refract_calib_file = os.path.join(raw_data_dir, 'refraction_FARO_times.csv')
            ocean_z_refract_calib_file = os.path.join(raw_data_dir, 'refraction_FARO_times.csv')
            # For now, use the dry scans to provide the shape of the final array to the processing steps, but do not leave it this way. eEventually you will need to the ZdataP file to construct refraction corrections.
            # postScanFile = os.path.join(matlab_data_dir, 'ZdataP.mat')
            postScanFile = dScanFile
            run_log_file = os.path.join(raw_data_dir, 'TDB_15_2_Log')
            # generates:
            # ocean levels for every hour
            # ocean levels during calibration scans
            # corrections for refraction.
            # saves them to disk
            topo_processing.generate_products(postScanFile, ocean_z_refract_calib_file, run_log_file, python_data_dir, _experiment)
            # then load them from disk because the function returns nothing
            _ocean_levels = np.load(ocean_level_file)
            _ref_corrections = np.load(refraction_correction_file)
        else:
            # otherwise, load them from disk
            _ocean_levels = np.load(ocean_level_file)
            _ref_corrections = np.load(refraction_correction_file)
            # raise Exception('called TDB_15_2 lookup function.')
        return _ocean_levels, _ref_corrections, _spurious_scans_dry, _spurious_scans_wet
    # This function is for TDB_12_1
    def tdb12_gen():
        print('generating necessary products...')
        _spurious_scans_dry = [328,474,656]
        _spurious_scans_wet = None        # print('called TDB_12_1 regeneration function!')
        if _regenerate:
            ocean_z_refract_calib_file = os.path.join(raw_data_dir, 'refraction_FARO_times.csv')
            postScanFile = dScanFile
            run_log_file = os.path.join(raw_data_dir, 'TDB_12_1 Ocean')
            topo_processing.generate_products(postScanFile, ocean_z_refract_calib_file, run_log_file, python_data_dir, _experiment)
            # then load them from disk because the function returns nothing
            _ocean_levels = np.load(ocean_level_file)
            _ref_corrections = np.load(refraction_correction_file)
            # raise Exception('called TDB_12_1 regeneration function!')
        else:
            # otherwise, load them from disk
            _ocean_levels = np.load(ocean_level_file)
            _ref_corrections = np.load(refraction_correction_file)
            # raise Exception('called TDB_12_1 lookup function.'
        return _ocean_levels, _ref_corrections, _spurious_scans_dry, _spurious_scans_wet

    def getAuxData(argument):
        experiments = {
            'tdb19': tdb19_gen,
            'tdb15': tdb15_gen,
            'tdb12': tdb12_gen
        }
        regenerateFunction = experiments.get(argument, 'nothing')
        return regenerateFunction()

    ocean_levels, ref_corrections, bad_dry_scans, bad_wet_scans = getAuxData(_experiment)
    _dims = dryScan.shape

    # print((len(bad_dry_scans), _dims[1],_dims[2]))
    # print((len(bad_wet_scans), _dims[1],_dims[2]))
    # print(bad_wet_scans)
    # # print(wetScan.shape)
    #
    # return dryScan.shape

    # if there is a scan where the data is suspect, replace it with zeros (275 in tdb19)
    if bad_dry_scans is not None:
        print('zeroing out bad dry scans')
        dryScan[bad_dry_scans,:,:] = np.zeros((len(bad_dry_scans), _dims[1], _dims[2]))
    if bad_wet_scans is not None:
        print('zeroing out bad wet scans')
        wetScan[bad_wet_scans,:,:] = np.zeros((len(bad_wet_scans), _dims[1], _dims[2]))

    # replace all scans that are missing (all zero) with the corresponding wet scan. if neither are ok, find nearest ok scans and find average.

    dryScan = topo_processing.swap_spare_maps(dryScan, wetScan)

    # This mask will be used a lot. The walls and the area behind the inlet are FALSE, the rest of the delta is TRUE.

    mask1 = topo_mask_gen.make_topo_mask1(dryScan.shape, _experiment)

    # make an empty array to hold all the scans once image processing has been finished.

    scans_cleaned = np.zeros(dryScan.shape)
    topset_masks = np.zeros(dryScan.shape, dtype = bool)
    dry_masks = np.zeros(dryScan.shape, dtype = bool)

    # return dryScan, mask1, ocean_levels, ref_corrections

    def tdb19_fill(_scan, _mask, _ocnz):
        return topo_processing.fill_topo(_scan, _mask, _ocnz)

    def tdb15_fill(_scan, _mask, _ocnz):
        return topo_processing.fill_topo(_scan, _mask, _ocnz)

    def tdb12_fill(_scan, _mask, _ocnz):
        return topo_processing.fill_topo(_scan, _mask, _ocnz, pit_thresh = 0.01, view_depth = 0, pitsize = 5)

    def getFillFunction(argument, _scan, _mask, _ocnz):
        experiments = {
            'tdb19': tdb19_fill,
            'tdb15': tdb15_fill,
            'tdb12': tdb12_fill
        }
        fillFunc = experiments.get(argument, 'nothing')
        return fillFunc(_scan, _mask, _ocnz)

    # main loop processing lidar images.
    print('cleaning lidar data...')
    for i in tqdm(np.arange(dryScan.shape[0])):

        scan = dryScan[i, :, :]
        scan[~mask1] = 0       # apply mask
        scan[scan > 0.65] = 0  # remove clearly spurious points

        # need data for ocean level
        ocean_z = ocean_levels[i]

        # fill the holes

        scan_filled, topset_mask = getFillFunction(_experiment, scan, mask1, ocean_z)

        dry_mask = topo_processing.get_dry_mask(scan, mask1, ocean_z)

        # correct for refraction
        # scan_corrected = topo_processing.correct_refraction(scan_filled, ocean_z, ref_corrections)

        # remove_mats and islands, smooth seabed
        # scan_mat_removal = topo_processing.remove_mats(scan_corrected, ocean_z, mask1)
        # scan_no_mats = scan_mat_removal[0]
        # mats_only = scan_mat_removal[1]
        dry_masks[i, :, :] = dry_mask
        scans_cleaned[i, :, :] = scan_filled
        topset_masks[i, :, :] = topset_mask

    # some ocean bottom values are still extra spiky and strange.
    # do additional processing to substitute zero values with previous values.

    # print('filling in ocean values...')
    # for i in tqdm(np.arange(1,scans_cleaned.shape[2])):
    #     scan1 = scans_cleaned[:,:,(i-1)]
    #     scan = scans_cleaned[:,:,i]
    #     zer = scan == 0
    #     scan[zer] = scan1[zer]

    # save cleaned-up data and masks to disk
    clean_scans_file = os.path.join(python_data_dir, 'processed_scans.npy')
    topset_masks_file = os.path.join(python_data_dir, 'topset_masks.npy')
    dry_masks_file = os.path.join(python_data_dir, 'dry_masks.npy')

    np.save(clean_scans_file, scans_cleaned)
    np.save(topset_masks_file, topset_masks)
    np.save(dry_masks_file, dry_masks)

    # scans_cleaned[scans_cleaned == 0] = np.nan

    synstrat = np.zeros(scans_cleaned.shape)

    tmax = scans_cleaned.shape[0]

    synstrat[tmax - 1, :, :] = scans_cleaned[tmax - 1, :, :]

    # process into synthetic stratigraphy
    print('generating synthetic stratigraphy')
    for i in tqdm(np.flip(np.arange(tmax-1))):
        surf1 = synstrat[i+1,:,:]
        surf2 = scans_cleaned[i,:,:]
        strat = np.fmin(surf1, surf2)
        synstrat[i,:,:] = strat

    # mask out portion of synstrat at every timestep which is not topset.

    synstrat[~topset_masks] = np.nan

    # save synthetic stratigraphy to disk
    datestring = datetime.today().strftime('%Y%m%d')
    synstrat_file = os.path.join(python_data_dir, 'synstrat_' + _experiment + '_' + datestring + '.npy')
    np.save(synstrat_file, synstrat)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Options for topo data processing')

    parser.add_argument('--root_dir', '-r', help = 'the path to the root directory of the project.')
    parser.add_argument('--generate', '-g', type = str2bool, default = False, help = 'generate intermediate data products from raw data.')
    # parser.add_argument('--spurious_scans', help = 'give as a tuple (x,x,...) all the scans that should be negated to zero.')
    parser.add_argument('--experiment', '-e', choices = ['tdb12', 'tdb15', 'tdb19'], help = 'a tag for the experiment id that is being processed, one of ("tdb19", "tdb15", or "tdb12")')

    args = parser.parse_args()
    # print(args.generate)

    process_lidar(args.root_dir, args.experiment, args.generate)
