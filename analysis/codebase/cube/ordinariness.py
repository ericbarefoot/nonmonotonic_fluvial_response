# functions to ordinariness in stratigraphic slices
# EA Barefoot
# Oct 2020

import numpy as np
from tqdm import tqdm


def detrend(etas):
    ''' takes eta as a T x X array
    '''
    pass


def get_original_percentiles(eta):
    ''' takes eta as a T x X array
    '''
    centeredEtas = eta - np.mean(eta, axis = 1)[:, np.newaxis]

    OGpercentiles = (centeredEtas.argsort(axis = 1).argsort(axis = 1) + 1) / lZ

    return OGpercentiles


def get_preserved_elevations(eta, synstrat):

    return eta == synstrat


def ordinariness(percentiles, preservedBool):

    preservedPercentiles = percentiles[preservedBool]

    ordinary = 1 - 2 * np.median(preservedPercentiles)
