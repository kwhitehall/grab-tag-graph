import sys

sys.path.insert(0, '/home/caocampb/PycharmProjects/grab-tag-graph/code/')

import iomethods
import variables

def test_read_MERG_and_write_MERG():
    '''
        Purpose:: Tests the functions 'read_MERG_pixel_file' and 'write_MERG_pixel_to_ncdf' from iomethods.py as a whole
        Notes:: Change the sys.path.insert(...) argument above to the path where iomethods.py is locateded
                (because iomethods_module_test.py is in a different folder than iomethods.py)
    '''
    user = variables.UserVariables(useJSON=False)

    lon, lat, temperatures = iomethods.read_MERG_pixel_file('/home/caocampb/PycharmProjects/grab-tag-graph/datadir/MERG/merg_2006091100_4km-pixel')

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

    iomethods.write_MERG_pixel_to_ncdf(lonDict, latDict, timeDict, ch4Dict, 'merg_2006091100_4km-pixel', user.DIRS['CEoriDirName'], globalAttrDict,
                             dimensionsDict)

if __name__ == 'main':
    test_read_MERG_and_write_MERG()


