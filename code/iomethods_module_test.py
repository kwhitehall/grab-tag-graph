import iomethods.py

if __name__ == '__main__':  # Testing for write_np_array_to_ncdf

    user = variables.UserVariables(useJSON=False)

    lon, lat, temperatures = read_MERG_pixel_file('/home/caocampb/PycharmProjects/grab-tag-graph/datadir/MERG/merg_2006091100_4km-pixel')

    lonDict = {"name": "longitude", "dataType": "double", "dimensions": ("longitude",), "units": "degrees_east", "long_name": "Longitude", "values": lon}

    latDict = {"name": "latitude", "dataType": "double", "dimensions": ("latitude",), "units": "degrees_north", "long_name": "Latitude", "values": lat}

    timeDict = {"name": "time", "dataType": "float", "dimensions": ("time",)}

    ch4Dict = {"name": "ch4", "dataType": "float", "dimensions": ("time", "latitude", "longitude"), "long_name": "IR BT (add 75 to this value)",
                                                                                                    "time_statistic": "instantaneous",
                                                                                                    "missing_value": float(330)}

    globalAttrDict = {"Conventions": "COARDS", "calendar": "standard", "comments": "File", "model": "geos/das",
                      "center": "gsfc"}

    dimensionsDict = {"time": None, "longitude": 9896, "latitude": 3298}

    write_MERG_pixel_to_ncdf(lonDict, latDict, timeDict, ch4Dict, 'mergFile', user.DIRS['CEoriDirName'], globalAttrDict,
                             dimensionsDict)
