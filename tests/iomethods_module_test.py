import sys
import timeit

from datetime import datetime

sys.path.insert(0, '../code/')

import iomethods
import variables
import utils
import numpy as np

def test_read_merg_and_write_merg():
    '''
        Purpose:: Tests the functions '_read_merg_file' and '_write_MERG_pixel_to_ncdf' from iomethods.py as a whole
        Notes:: Change the sys.path.insert(...) argument above to the path where iomethods.py is located
                (because iomethods_module_test.py is in a different folder than iomethods.py)
    '''

    user = variables.UserVariables(useJSON=False)

    lon, lat, temperatures = iomethods._read_merg_file('../datadir/MERG/merg_2006091100_4km-pixel', shape=(2, 3298, 9896), offset=75.)

    # Generate lon and lat coordinates
    lon = np.arange(0.0182, 360., 0.036378335, dtype=np.float)
    lat = np.arange(59.982, -60., -0.036383683, dtype=np.float)

    lonDict = {"name": "longitude", "dataType": "double", "dimensions": ("longitude",), "units": "degrees_east", "long_name": "Longitude", "values": lon}

    latDict = {"name": "latitude", "dataType": "double", "dimensions": ("latitude",), "units": "degrees_north", "long_name": "Latitude", "values": lat}

    timeDict = {"name": "time", "dataType": "float", "dimensions": ("time",)}

    ch4Dict = {"name": "ch4", "dataType": "float", "dimensions": ("time", "latitude", "longitude"), "long_name": "IR BT (add 75 to this value)",
                                                                                                    "time_statistic": "instantaneous",
                                                                                                    "missing_value": float(330),
                                                                                                    "values": temperatures}

    globalAttrDict = {"Conventions": "COARDS", "calendar": "standard", "comments": "File", "model": "geos/das",
                      "center": "gsfc"}

    dimensionsDict = {"time": None, "longitude": 9896, "latitude": 3298}

    iomethods.write_merg_to_ncdf(lonDict, latDict, timeDict, ch4Dict, 'merg_2006091100_4km-pixel', user.DIRS['CEoriDirName'], globalAttrDict,
                             dimensionsDict)

def test_read_netCDF_to_array():
    # Test for TRMM file, extracting the 'irp' variable
    minDate = datetime(2005, 1, 1)
    maxDate = datetime(2007, 1, 1)

    trimmedData, times, trimmedLats, trimmedLons = iomethods.read_netCDF_to_array('../datadir/TRMM/3B42.20090831.00.7A.nc',
                                                                           'trmm', 'irp', minDate, maxDate, 10, 15, 10, 15)

    # minDate = datetime(2009, 8, 21)
    # maxDate = datetime(2009, 8, 31)
    #
    # trimmedData, times, trimmedLats, trimmedLons = read_netCDF_to_array(
    #         '/Users/diegovonbeck/grab-tag-graph-bak/datadir/TRMM/3B42.20090831.00.7A.nc',
    #          'trmm', 'irp', minDate, maxDate, 100, 15, 10, 15)

    # Compare these with known values
    print trimmedData
    print times
    print trimmedLats
    print trimmedLons

def performanceTestDataSubset1():
    lats = np.arange(0, 60., 0.05, dtype=np.float)

    trimmedLatsIndices = [i for i in range(len(lats))
                          if lats[i] <= 50 and lats[i] >= 30]

    trimmedLatsStart = trimmedLatsIndices[0]
    trimmedLatsEnd = trimmedLatsIndices[-1]

    trimmedLats = lats[trimmedLatsEnd:trimmedLatsStart]


def performanceTestDataSubset2():
    lats = np.arange(0, 60., 0.05, dtype=np.float)

    latminNETCDF = utils.find_nearest(lats, float(30))
    latmaxNETCDF = utils.find_nearest(lats, float(50))

    latminIndex = (np.where(lats == latminNETCDF))[0][0]
    latmaxIndex = (np.where(lats == latmaxNETCDF))[0][0]

    trimmedLats = lats[latminIndex:latmaxIndex]


if __name__ == '__main__':
    
    start_time = timeit.default_timer()
    performanceTestDataSubset1()
    print(timeit.default_timer() - start_time)

    start_time = timeit.default_timer()
    performanceTestDataSubset2()
    print(timeit.default_timer() - start_time)

    test_read_netCDF_to_array ()


