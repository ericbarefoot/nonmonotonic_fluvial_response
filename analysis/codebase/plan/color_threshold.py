# Function that takes overhead images of TDB experiments, applies an appropriate mask,
# then finds a guess for the topset via image processing, and saves the mask as well
# as the ndwi for thresholding and secondary processing.

# Eric Barefoot
# Nov 2019

# import the necessary packages
import os
import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import rc
# rc('image', cmap="Greys_r")
from matplotlib.path import Path
from scipy import ndimage
from skimage.io import imread, imsave
from skimage import morphology, filters, data, transform, img_as_bool, measure, filters, img_as_ubyte, img_as_float
from scipy.ndimage import rotate
from skimage import transform as trans
from tqdm import tqdm

def mask_maker(mask, experiment):

    # make a mask to clip out the area not around the delta.

    # These parameters are for TDB_19_1
    def tdb19_params():
        W = 744
        H = int(W * 4/3)
        widthBoundMult = 0.32
        heightBoundMult = 0.175
        rotation = 135
        scale = 1.6
        flip = False
        return rotation, widthBoundMult, heightBoundMult, scale, H, W, flip

    # These parameters are for TDB_15_2
    # will also have to deal with the fact that theres .json files.
    def tdb15_params():
        # raise Exception('called function to make parameters for TDB 15')
        W = 732
        H = int(W * 4/3)
        widthBoundMult = 0.03
        heightBoundMult = 0.14
        rotation = -44
        scale = 1.6
        flip = False
        return rotation, widthBoundMult, heightBoundMult, scale, H, W, flip

    # These parameters are for TDB_12_1
    def tdb12_params():
        # raise Exception('called function to make parameters for TDB 12')
        W = 600
        H = W
        widthBoundMult = 0.02
        heightBoundMult = 0.08
        rotation = -2
        scale = 0.88
        flip = False
        return rotation, widthBoundMult, heightBoundMult, scale, H, W, flip

    experiments = {
        'tdb19': tdb19_params,
        'tdb15': tdb15_params,
        'tdb12': tdb12_params
    }

    def get_parameters(argument):
        params = experiments.get(argument, 'nothing')
        return params()

    # return get_parameters(experiment)
    angle, wMult, hMult, scaleFactor, newH, newW, flipBool = get_parameters(experiment)

    if experiment == 'tdb12':

        firstpad = int((mask.shape[0] - mask.shape[1]) / 2)
        padone = np.pad(mask, ((0, 0), (firstpad, firstpad)), 'constant')
        paddMask = trans.resize(padone, (newW, newW), order = 0)

        flipMask = rotate(np.flip(paddMask, 1), 90, order = 0)
        scaleMask = trans.rescale(flipMask, scaleFactor, order = 0, multichannel = False)

        mw, mh = scaleMask.shape

        flipScale = rotate(scaleMask, -2, order = 0)

        secondpad1 = int(mw * wMult)
        secondpad2 = int(mh * hMult)
        padtwo = np.pad(flipScale, ((secondpad1, 0), (secondpad2, 0)), 'constant')
        padthree = np.pad(padtwo, ((0, newW - padtwo.shape[0]), (0, newH - padtwo.shape[1])), 'constant')
        topset_mask = padthree.astype('bool')
    else:
        rot_mask = rotate(mask, angle, order = 0)

        mask_sized = trans.resize(rot_mask, (newW, newW), order = 0)

        padd = int((newH - newW) / 2)

        paddMask = np.pad(mask_sized, ((0, 0), (padd, padd)), 'constant')

        upscaleMask = trans.rescale(paddMask, scaleFactor, order = 0, multichannel = False)

        mw, mh = upscaleMask.shape

        lboundW = int(mw * wMult)
        uboundW = lboundW + newW
        lboundH = int(mh * hMult)
        uboundH = lboundH + newH

        if flipBool:
            flipMask = np.flip(upscaleMask, 1)
        else:
            flipMask = upscaleMask

        topset_mask = flipMask[lboundW:uboundW, lboundH:uboundH].astype('bool')

    return topset_mask

# def rescale_images(images):
#
#     experiments = {
#         'tdb19': (744, 992),
#         'tdb15': (732, 976),
#         'tdb12': (600, 600)
#     }
#
#     def get_parameters(argument):
#         params = experiments.get(argument, 'nothing')
#         return params
#
#     newW, newH = get_parameters(experiment)
#
#     stacked_images = np.zeros((newW, newH, 3, images.shape[-1]))
#
#     imgRS = trans.resize(img, (newW, newH))
#     stacked_images[:,:,:,ind] = imgRS
#
#     return stacked_images

def get_ndwi(image, mask):

    # split the images into bands

    image[~mask] = np.nan

    r = image[...,0]
    g = image[...,1]
    b = image[...,2]

    # r[~mask] = np.nan
    # g[~mask] = np.nan
    # b[~mask] = np.nan

    # make a water blueness index to threshold the image.

    np.seterr(divide='ignore', invalid='ignore')

    check = np.logical_or(g > 0, r > 0)

    ndwi = np.where(check, (g - r) / (b * (g + r)), np.nan)

    return ndwi

def isolate_channels(ndwi, _experiment, dilate_star = 1):

    experiments = {
        'tdb19': (1.5, 1000, 100),
        'tdb15': (0.75, 1000, 100),
        'tdb12': (1.5, 1000, 100)
    }

    def get_threshs(argument):
        return experiments.get(argument, 'nothing')

    ndwiThresh, smObj, smHole = get_threshs(_experiment)

    chan = ndwi > (ndwiThresh * np.nanmean(ndwi))

    channel_mask = np.zeros(chan.shape, dtype = 'bool')

    chan_rem = morphology.remove_small_objects(chan, min_size = smObj)
    chan_fill = morphology.remove_small_holes(chan_rem, area_threshold = smHole)
    # channel_mask = morphology.erosion(chan_fill, selem = morphology.star(a = dilate_star))

    return chan_fill
