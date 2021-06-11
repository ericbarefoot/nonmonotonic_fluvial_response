from matplotlib.path import Path
import numpy as np

def make_topo_mask1(_dims, _experiment, _data_map = None):

    height = _dims[1]
    width = _dims[2]

    def tdb19_mask():
        polygon1 = [(29,27),(74,25),(75,66),(65,77),(25,77)]
        polygon1 = [(438,401),(378,425),(339,474),(282,550),(376,550)]
        polygon1 = [(0,0),(162,0),(159,7),(139,14),(119,36),(98,58),(81,75),(57,98),(32,123),(11,148),(0,150)]
        poly1_path = Path(polygon1)

        x, y = np.mgrid[:height, :width]
        coors=np.hstack((x.reshape(-1, 1), y.reshape(-1,1)))

        _mask = ~poly1_path.contains_points(coors).reshape(height, width)
        _mask[:,0:6] = 0
        _mask[0:6,:] = 0
        return _mask

    def tdb15_mask():
        polygon1 = [
        (0,0),
        (162,0),
        (159,7),
        (139,14),
        (119,36),
        (98,58),
        (81,75),
        (57,98),
        (32,123),
        (11,148),
        (0,150)]
        polygon2 = [
        (400,0),
        (400,100),
        (415,180),
        (425,300),
        (280,440),
        (145,400),
        (0,350),
        (0, width),
        (height, width),
        (height, 0)]

        poly1_path = Path(polygon1)
        poly2_path = Path(polygon2)

        x, y = np.mgrid[:height, :width]
        coors = np.hstack((x.reshape(-1, 1), y.reshape(-1,1)))

        _mask1 = ~poly1_path.contains_points(coors).reshape(height, width)
        _mask2 = ~poly2_path.contains_points(coors).reshape(height, width)
        _mask = np.logical_and(_mask1, _mask2)
        _mask[:,0:7] = 0
        _mask[0:7,:] = 0
        return _mask

    def tdb12_mask():
        polygon1 = [
        (160,0),
        (86,77),
        (0,165),
        (0,380),
        (86,422),
        (321,337),
        (width,210),
        (width,0)]

        poly1_path = Path(polygon1)
        x, y = np.mgrid[:height, :width]
        coors = np.hstack((x.reshape(-1, 1), y.reshape(-1,1)))

        _mask = poly1_path.contains_points(coors).reshape(height, width)
        _mask[:,0:5] = 0
        _mask[0:5,:] = 0
        return _mask

    def getMask(argument):
        experiments = {
            'tdb19': tdb19_mask,
            'tdb15': tdb15_mask,
            'tdb12': tdb12_mask
        }
        getMaskFun = experiments.get(argument, 'nothing')
        return getMaskFun()

    mask1 = getMask(_experiment)

    return mask1
