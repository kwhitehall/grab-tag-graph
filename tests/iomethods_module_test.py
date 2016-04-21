import sys
import timeit

from datetime import datetime

sys.path.insert(0, '../code/')

import iomethods
import variables
import utils
import unittest
import numpy as np

class TestIO(unittest.TestCase):
    def test_read_netCDF_to_array(self):
        minDate = datetime(2009, 8, 31)
        maxDate = datetime(2009, 8, 31)

        trimmedData, times, trimmedLats, trimmedLons = iomethods.read_netCDF_to_array('../datadir/TRMM/3B42.20090831.00.7A.nc',
                                                                               'irp', minDate, maxDate, 10, 15, 10, 15)

        #TODO, validate trimmedData
        '''
        print trimmedData
        expectedTrimmedData = [[[0.17999999, 0.44, 0.57999998, 0.57999998, 0.39999998, 0.13, 0., 0.25999999, 1.81999993, 3.36999989, 4.9000001, 5.42999983, 3.86999989, 1.42999995, 0.44999999, 0.09, 0., 0.32999998, 0., 0., 0.],
                                [0.11, 0.19999999, 0.22999999, 0.29999998, 0.13, 0., 0., 0.17999999, 1.14999998, 2.54999995, 4.00999975, 4.00999975, 3.12999988, 0.81999999, 0., 0., 0.09, 0., 0., 0., 0.],
                                [0., 0.14999999, 0.11, 0., 0., 0., 0., 0.19999999, 0.53999996, 1.12, 1.62, 2.22000003, 1.31999993, 0.39999998, 0., 0., 0.25, 0., 0., 0., 0.],
                                [0.61000001, 0.29999998, 0., 0., 0., 0., 0., 0.25, 0.22999999, 0., 0., 0.34, 0.23999999, 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0.19, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]]]
        '''
        expectedTimes = [datetime(2009, 8, 31, 0, 0)]
        expectedLats = [9.875, 10.125, 10.375, 10.625, 10.875, 11.125, 11.375, 11.625, 11.875, 12.125, 12.375,\
                        12.625, 12.875, 13.125, 13.375, 13.625, 13.875, 14.125, 14.375, 14.625, 14.875]
        expectedLons = [9.875, 10.125, 10.375, 10.625, 10.875, 11.125, 11.375, 11.625, 11.875, 12.125, 12.375,\
                        12.625, 12.875, 13.125, 13.375, 13.625, 13.875, 14.125, 14.375, 14.625, 14.875]

        #self.assertTrue((expectedTrimmedData == trimmedData).all(), "Trimmed and expected data differ!")
        self.assertEqual(expectedTimes, times, "Returned time and expected time differ!")
        self.assertTrue((expectedLats == trimmedLats).all(), "Trimmed and expected Lat ranges differ!")
        self.assertTrue((expectedLons == trimmedLons).all(), "Trimmed and expected Lon ranges differ!")

    def test_read_merg_and_write_merg(self):
        '''
            Purpose:: Tests the functions '_read_merg_file' and '_write_MERG_pixel_to_ncdf' from iomethods.py as a whole
        '''
        FILE = "../datadir/MERG/merg_2006091100_4km-pixel"
        user = variables.UserVariables(useJSON=False)
        self.assertIsNotNone(user, "Invalid user variable given")

        #lon, lat, temperatures = iomethods._read_merg_file('../datadir/MERG/merg_2006091100_4km-pixel', shape=(2, 3298, 9896), offset=75.)
        data = iomethods._read_merg_file(FILE, shape=(2, 3298, 9896), offset=75.)

        # Generate lon and lat coordinates
        lon = np.arange(0.0182, 360., 0.036378335, dtype=np.float)
        lat = np.arange(59.982, -60., -0.036383683, dtype=np.float)

        # Generate dictionaries with various dimensions, attributes, and variables
        lonDict = {"name": "longitude", "dataType": "double", "dimensions": ("longitude",), "units": "degrees_east", "long_name": "Longitude", "values": lon}
        latDict = {"name": "latitude", "dataType": "double", "dimensions": ("latitude",), "units": "degrees_north", "long_name": "Latitude", "values": lat}
        timeDict = {"name": "time", "dataType": "float", "dimensions": ("time",)}
        ch4Dict = {"name": "ch4", "dataType": "float", "dimensions": ("time", "latitude", "longitude"), "long_name": "IR BT (add 75 to this value)",
                                                                                                        "time_statistic": "instantaneous",
                                                                                                        "missing_value": float(330),
                                                                                                        "values": data}
        globalAttrDict = {"Conventions": "COARDS", "calendar": "standard", "comments": "File", "model": "geos/das",
                          "center": "gsfc"}
        dimensionsDict = {"time": None, "longitude": 9896, "latitude": 3298}

        output = iomethods.write_merg_to_ncdf(lonDict, latDict, timeDict, ch4Dict, 'merg_2006091100_4km-pixel', user.DIRS['CEoriDirName'],
                                              globalAttrDict, dimensionsDict)
        self.assertIsNotNone(output, "Returned output is empty!")
        
        #TODO validate contents of file, not that it just exists
        self.assertEqual(output, FILE + ".nc", "Output file name differs than expected!")


def list_comprehension_trimmer():
    '''
    Purpose:: Trim's array of numbers outside a certain range by using a list comprehension
              It's to be compared with find_nearest_trimmer which does the same thing in a different way

    '''
    lats = np.arange(0, 60., 0.05, dtype=np.float)

    trimmedLatsIndices = [i for i in range(len(lats))
                          if lats[i] <= 50 and lats[i] >= 30]

    trimmedLatsStart = trimmedLatsIndices[0]
    trimmedLatsEnd = trimmedLatsIndices[-1]

    trimmedLats = lats[trimmedLatsEnd:trimmedLatsStart]


def find_nearest_trimmer():
    '''
    Purpose:: Trim's array of numbers outside a certain range by using utils.find_nearest() and np.where()
              It's to be compared with list_comprehension_trimmer which does the same thing in a different way

    '''
    lats = np.arange(0, 60., 0.05, dtype=np.float)

    latminNETCDF = utils.find_nearest(lats, float(30))
    latmaxNETCDF = utils.find_nearest(lats, float(50))

    latminIndex = (np.where(lats == latminNETCDF))[0][0]
    latmaxIndex = (np.where(lats == latmaxNETCDF))[0][0]

    trimmedLats = lats[latminIndex:latmaxIndex]


if __name__ == '__main__':
    '''
    start_time = timeit.default_timer()
    performanceTestDataSubset1()
    print(timeit.default_timer() - start_time)

    start_time = timeit.default_timer()
    performanceTestDataSubset2()
    print(timeit.default_timer() - start_time)
    '''
    unittest.main()


