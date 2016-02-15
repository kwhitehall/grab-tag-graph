import glob
import itertools
import json
import os
import time
from datetime import timedelta
from multiprocessing import Manager
from multiprocessing import Pool

import networkx as nx
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset, date2num
from scipy import ndimage

# GTG modules
import utils

NUM_IMAGE_WORKERS = 2  # Number of workers to send off for extracting CE in the independent image frames
P_TIME = 0

# ------------------------ End GLOBAL VARS -------------------------
# ************************ Begin Functions *************************
# **********************************************************************************************************************
manager = Manager()
varsDict = manager.dict()


# Callable object that is passed to the Pool map (so that the find_cloud_elements function can be called many times in parallel) - Added by Gabriel Mel
class CeFinder(object):
    '''

    '''
    def __init__(self, timelist, userVariables, TRMMdirName=None):
        self.timelist = timelist
        self.userVariables = userVariables
        self.TRMMdirName = TRMMdirName

    def __call__(self, t):
        return find_single_frame_cloud_elements(t, varsDict['images'], self.timelist,
                                                varsDict['lat'], varsDict['lon'], self.userVariables, self.TRMMdirName)
# **********************************************************************************************************************


def find_cloud_elements(mergImgs, timelist, lat, lon, userVariables, graphVariables, TRMMdirName=None):
    '''
    Purpose: parallelizes the process to determine the contiguous boxes for a given time in the satellite images 
             i.e. each frame using scipy ndimage package

    Inputs:: mergImgs: masked numpy array in (time,lat,lon),T_bb representing the satellite data. This is masked based 
        on the maximum acceptable temperature, T_BB_MAX
        timelist: a list of python datatimes
        lat: a 2D array of the set of latitudes from the files opened
        lon: a 2D array of the set of longitudes from the files opened
        userVariables:
        graphVariables:
        TRMMdirName (optional): string representing the path where to find the TRMM datafiles

    Returns::
        a Networkx directed graph where each node contains the information in cloudElementDict
        The nodes are determined according to the area of contiguous squares. The nodes are linked through weighted
        edges.

    Assumptions::
        Assumes we are dealing with MERG data which is 4kmx4km resolved, thus the smallest value
        required according to Vila et al. (2008) is 2400km^2
        therefore, 2400/16 = 150 contiguous squares
    '''

    varsDict['images'] = mergImgs
    varsDict['lat'] = lat
    varsDict['lon'] = lon

    global LAT
    LAT = lat
    global LON
    LON = lon

    p = Pool(NUM_IMAGE_WORKERS)
    image_proc_start = time.time()

    results = p.map(CeFinder(timelist, userVariables, TRMMdirName), xrange(mergImgs.shape[0]))

    return assemble_graph(results, userVariables, graphVariables)
# **********************************************************************************************************************


def assemble_graph(results, userVariables, graphVariables):
    '''
    Purpose:: to build the full graph from all the data being considered. 
        i.e. Once the images have been processed in a parallelized manner, put the results together 

    Inputs:: results: a list of two items
               (1) a list of dictionaries for all cloud elements found
               (2) a string with the information for all cloud elements meeting the size and T criteria 
             userVariables: 
             graphVariables: 

    Outputs:: CLOUD_ELEMENT_GRAPH: a Networkx directed graph where each node contains the information in cloudElementDict
        The nodes are determined according to the area of contiguous squares. The nodes are linked through weighted
        edges.
        ceCountTuple: a tuple of list of ints representing (totalCEsList, acceptedCEsList) where 
            totalCEsList is a list of the sum of all CEs identified for each time between starttime and endtime and 
            acceptedCEsList is a list of the sum of CEs connected to another frame for each time between starttime and endtime in the graph
    
    Assumptions::
    '''

    totalCEs = 0
    totalCEsList = []
    acceptedCEs = 0
    acceptedCEsList = []

    edgeWeight = [1, 2, 3]  # weights for the graph edges
    cloudElementsJSON = []  # list of the key, value objects associated with a CE in the graph
    seenNode = []
    filenameJSON = userVariables.DIRS['mainDirStr'] + '/textFiles/graphJSON.txt'

    # openfile for storing ALL cloudElement information
    cloudElementsFile = open((userVariables.DIRS['mainDirStr'] + '/textFiles/cloudElements.txt'), 'wb')
    # openfile for storing cloudElement information meeting user criteria i.e. MCCs in this case
    cloudElementsUserFile = open((userVariables.DIRS['mainDirStr'] + '/textFiles/cloudElementsUserFile.txt'), 'w')

    for ce in results[0][0]:
        if ce['uniqueID'] not in dict(enumerate(zip(*seenNode))):
            totalCEs += 1
            graphVariables.CLOUD_ELEMENT_GRAPH.add_node(ce['uniqueID'], ce)
            seenNode.append(ce['uniqueID'])
            cloudElementsUserFile.write('Time is: %s' % (str(ce['cloudElementTime'])))
            cloudElementsUserFile.write('\nceUniqueID is: %s' % ce['uniqueID'])
            cloudElementsUserFile.write('\nCenter (lat,lon) is: %s' % ce['cloudElementCenter'])
            cloudElementsUserFile.write('\nArea is: %.4f km^2' % ce['cloudElementArea'])
            cloudElementsUserFile.write('\nAverage brightness temperature is: %.4f' % ce['cloudElementBTavg'])
            cloudElementsUserFile.write('\nMin brightness temperature is: %.4f K' % ce['cloudElementTmin'])
            cloudElementsUserFile.write('\nMax brightness temperature is: %.4f K' % ce['cloudElementTmax'])
            cloudElementsUserFile.write('\nBrightness temperature variance is: %.4f K' % ce['cloudElementBTvar'])
            cloudElementsUserFile.write('\nConvective fraction is: %.4f ' % ce['cloudElementCF'])
            cloudElementsUserFile.write('\nEccentricity is: %.4f ' % ce['cloudElementEccentricity'])

            cloudElementsFile.write(results[0][1])

    for t in xrange(1, len(results)):
        currFrameCEs = results[t][0]
        prevFrameCEs = results[t-1][0]
        ceNum = 0

        for ce in currFrameCEs:
            if ce['uniqueID'] not in dict(enumerate(zip(*seenNode))):
                graphVariables.CLOUD_ELEMENT_GRAPH.add_node(ce['uniqueID'], ce)
                seenNode.append(ce['uniqueID'])
                edges = []
                cloudElementsUserFile.write('\n\nTime is: %s' % (str(ce['cloudElementTime'])))
                cloudElementsUserFile.write('\nceUniqueID is: %s' % ce['uniqueID'])
                cloudElementsUserFile.write('\nCenter (lat,lon) is: %s' % ce['cloudElementCenter'])
                cloudElementsUserFile.write('\nArea is: %.4f km^2' % ce['cloudElementArea'])
                cloudElementsUserFile.write('\nAverage brightness temperature is: %.4f' % ce['cloudElementBTavg'])
                cloudElementsUserFile.write('\nMin brightness temperature is: %.4f K' % ce['cloudElementTmin'])
                cloudElementsUserFile.write('\nMax brightness temperature is: %.4f K' % ce['cloudElementTmax'])
                cloudElementsUserFile.write('\nBrightness temperature variance is: %.4f K' % ce['cloudElementBTvar'])
                cloudElementsUserFile.write('\nConvective fraction is: %.4f ' % ce['cloudElementCF'])
                cloudElementsUserFile.write('\nEccentricity is: %.4f ' % ce['cloudElementEccentricity'])

                cloudElementsFile.write(results[t][1])
                totalCEs += 1

                for cloudElementDict in prevFrameCEs:
                    percentageOverlap, areaOverlap = cloud_element_overlap(ce['cloudElementLatLon'],
                                                                           cloudElementDict['cloudElementLatLon'], userVariables.XRES, userVariables.YRES)

                    if percentageOverlap >= 0.95:
                        graphVariables.CLOUD_ELEMENT_GRAPH.add_edge(cloudElementDict['uniqueID'], ce['uniqueID'], weight=graphVariables.edgeWeight[0])
                        edges.append(cloudElementDict['uniqueID'])
                        acceptedCEs += 1

                    elif percentageOverlap >= 0.90 and percentageOverlap < 0.95:
                        graphVariables.CLOUD_ELEMENT_GRAPH.add_edge(cloudElementDict['uniqueID'], ce['uniqueID'], weight=graphVariables.edgeWeight[1])
                        edges.append(cloudElementDict['uniqueID'])
                        acceptedCEs += 1

                    elif areaOverlap >= userVariables.MIN_OVERLAP:
                        graphVariables.CLOUD_ELEMENT_GRAPH.add_edge(cloudElementDict['uniqueID'], ce['uniqueID'], weight=graphVariables.edgeWeight[2])
                        edges.append(cloudElementDict['uniqueID'])
                        acceptedCEs += 1

                if edges:
                    cloudElementsJSON.append({'cloudElement': ce['uniqueID'], 'time': str(ce['cloudElementTime']),
                        'area': ce['cloudElementArea'], 'Tmax': ce['cloudElementTmax'], 'Tmin': ce['cloudElementTmin'], 'center': ce['cloudElementCenter'],
                        'shape': ce['cloudElementEccentricity'], 'cloudElementLatLonBox': ce['cloudElementLatLonBox'], 'convective_fraction': ce['cloudElementCF'], 'edges': edges})

        acceptedCEsList.append(acceptedCEs)
        acceptedCEs = 0
        totalCEsList.append(totalCEs)
        totalCEs = 0

    # Close info files
    cloudElementsFile.close()
    cloudElementsUserFile.close()

    # Write to JSON file
    with open(filenameJSON, 'w+') as f:
        json.dump(cloudElementsJSON, f)

    # Clean up graph - remove parent and childless nodes
    outAndInDeg = graphVariables.CLOUD_ELEMENT_GRAPH.degree_iter()
    toRemove = [node[0] for node in outAndInDeg if node[1] < 1]
    graphVariables.CLOUD_ELEMENT_GRAPH.remove_nodes_from(toRemove)

    return graphVariables.CLOUD_ELEMENT_GRAPH, (totalCEsList, acceptedCEsList)
# **********************************************************************************************************************


def find_single_frame_cloud_elements(t, mergImgs, timelist, lat, lon, userVariables, TRMMdirName=None):
    '''
    Purpose:: Determines the contiguous boxes for a given time of the satellite images i.e. each frame
        using scipy ndimage package in a parallelized manner. This function is called by many workers. 

    Inputs:: mergImgs: masked numpy array in (time,lat,lon),T_bb representing the satellite data. This is masked based 
        on the maximum acceptable temperature, T_BB_MAX
        timelist: a list of python datatimes
        lat: a 2D array of the set of latitudes from the files opened
        lon: a 2D array of the set of longitudes from the files opened
        TRMMdirName (optional): string representing the path where to find the TRMM datafiles

    Returns:: a list of two items
               (1) a list of dictionaries for all cloud elements found
               cloudElementDict = {'uniqueID': unique tag for this CE,
                            'cloudElementTime': time of the CE,
                            'cloudElementLatLon': (lat,lon,value) of MERG data of CE,
                            'cloudElementCenter':list of floating-point [lat,lon] representing the CE's center
                            'cloudElementArea':floating-point representing the area of the CE,
                            'cloudElementEccentricity': floating-point representing the shape of the CE,
                            'cloudElementTmax':integer representing the maximum Tb in CE,
                            'cloudElementTmin': integer representing the minimum Tb in CE,
                            'cloudElementCF': a floating-point representing the convective_fraction
                            'cloudElementLatLonBox': a list of floating-point representing the corners of the bounding box
                                                    of the CE [min_cloudElementLon, min_cloudElementLat, max_cloudElementLon,
                                                    max_cloudElementLat]
                            'cloudElementBTvar': floating-point representing the variance of BT across the CE
                            'cloudElementBTavg': floating-point representing the average of BT across the CE
                            'cloudElementPrecipTotal':floating-point representing the sum of all rainfall in CE if
                                                      TRMMdirName entered,
                            'cloudElementLatLonTRMM':(lat,lon,value) of TRMM data in CE if TRMMdirName entered,
                            'TRMMArea': floating-point representing the CE if TRMMdirName entered,
                            'CETRMMmax':floating-point representing the max rate in the CE if TRMMdirName entered,
                            'CETRMMmin':floating-point representing the min rate in the CE if TRMMdirName entered}
                            
               (2) a string with the information for all cloud elements meeting the size and T criteria
        
    Assumptions::
        Assumes we are dealing with MERG data which is 4kmx4km resolved, thus the smallest value
        required according to Vila et al. (2008) is 2400km^2
        therefore, 2400/16 = 150 contiguous squares
    '''
    single_frame_start = time.time()

    global LAT
    LAT = lat
    global LON
    LON = lon

    # global P_TIME

    cloudElementsJSON = []  # list of the key, value objects associated with a CE in the graph
    edges = []              # list of the nodes connected to a given CE
    latLonBox = []          # list of extreme points for the min fitting box around a CE [lon_min, lat_min, lon_max, lat_max]
    shape = 0               # max(num of non-zero boxes in lat, num of non-zero boxes in lon) for the CE
    cf = 0.0                # convective fraction
    tmin = 0.0              # IR Tmin for the CE
    tmax = 0.0              # IR Tmax for the CE
    varBT = 0.0             # IR variance
    avgBT = 0.0             # IR average

    # frame = ma.empty((1, mergImgs.shape[1], mergImgs.shape[2]))
    ceCounter = 0
    frameceCounter = 0
    frameNum = 0
    cloudElementEpsilon = 0.0
    cloudElementDict = {}
    cloudElementCenter = []     # list with two elements [lat,lon] for the center of a CE
    prevFrameCEs = []           # list for CEs in previous frame
    currFrameCEs = []           # list for CEs in current frame
    cloudElementLat = []        # list for a particular CE's lat values
    cloudElementLon = []        # list for a particular CE's lon values
    cloudElementLatLons = []    # list for a particular CE's (lat,lon) values
    allCloudElementDicts = []

    temporalRes = 3  # TRMM data is 3 hourly
    precipTotal = 0.0
    ceTRMMList = []
    precip = []

    maxCEprecipRate = 0.0
    minCEprecipRate = 0.0

    nygrd = len(LAT[:, 0])
    nxgrd = len(LON[0, :])

    cloudElementsFileString = ''

    # Determine contiguous locations with temperature below the warmest temp i.e. cloudElements in each frame
    frame, ceCounter = ndimage.measurements.label(mergImgs[t, :, :], structure=userVariables.STRUCTURING_ELEMENT)
    frameceCounter = 0
    frameNum = t + 1

    # If other dataset (TRMM) assumed to be a precipitation dataset was entered,
    # first do the regridding to avoid duplicating this effort
    # NB in the TRMM files the info is hours since the time thus 00Z file has in 01, 02 and 03 times
    if TRMMdirName:
        #  ------------------TRMM stuff -------------------------------------------------
        fileDate = ((str(timelist[t])).replace(' ', '')[:-8]).replace('-', '')
        fileHr1 = (str(timelist[t])).replace(' ', '')[-8:-6]

        if int(fileHr1) % temporalRes == 0:
            fileHr = fileHr1
        else:
            fileHr = (int(fileHr1) / temporalRes) * temporalRes
        if fileHr < 10:
            fileHr = '0'+str(fileHr)
        else:
            str(fileHr)

        # Open TRMM file for the resolution info and to create the appropriate sized grid
        TRMMfileName = TRMMdirName + '/3B42.' + fileDate + '.' + str(fileHr) + '.7A.nc'

        TRMMData = Dataset(TRMMfileName, 'r', format='NETCDF4')
        precipRate = TRMMData.variables['pcp'][:, :, :]
        latsrawTRMMData = TRMMData.variables['latitude'][:]
        lonsrawTRMMData = TRMMData.variables['longitude'][:]
        lonsrawTRMMData[lonsrawTRMMData > 180] = lonsrawTRMMData[lonsrawTRMMData > 180] - 360.
        LONTRMM, LATTRMM = np.meshgrid(lonsrawTRMMData, latsrawTRMMData)

        precipRateMasked = ma.masked_array(precipRate, mask=(precipRate < 0.0))
        # ---------regrid the TRMM data to the MERG dataset ----------------------------------
        # regrid using the do_regrid stuff from the Apache OCW
        regriddedTRMM = ma.zeros((0, nygrd, nxgrd))
        regriddedTRMM = utils.do_regrid(precipRateMasked[0, :, :], LATTRMM, LONTRMM, LAT, LON, order=1,
                                        mdi=-999999999)
        TRMMData.close()

    # For each of the areas identified in the IR data, check to determine valid CEs via an area and T requirement
    for count in xrange(ceCounter):
        # [0] is time dimension. Determine the actual values from the data
        # loc is a masked array
        try:
            loc = ndimage.find_objects(frame == (count+1))[0]

        except Exception, e:
            print 'Error is ', e
            continue

        cloudElement = mergImgs[t, :, :][loc]
        labels, _ = ndimage.label(cloudElement)

        # Determine the true lats and lons for this particular CE
        cloudElementLat = LAT[loc[0], 0]
        cloudElementLon = LON[0, loc[1]]

        # Determine number of boxes in this cloudelement
        numOfBoxes = np.count_nonzero(cloudElement)
        cloudElementArea = numOfBoxes * userVariables.XRES * userVariables.YRES

        # If the area is greater than the area required, or if the area is smaller than the suggested area,
        # check if it meets a convective fraction requirement consider as CE

        if cloudElementArea >= userVariables.AREA_MIN or (cloudElementArea < userVariables.AREA_MIN and
                    ((ndimage.minimum(cloudElement, labels=labels)) /
                    float((ndimage.maximum(cloudElement, labels=labels)))) < userVariables.CONVECTIVE_FRACTION):

            # Get some time information and labeling info
            frameceCounter += 1
            ceUniqueID = 'F' + str(frameNum) + 'CE' + str(frameceCounter)
            # -------------------------------------------------
            # Textfile name for accessing CE data using MATLAB code
            # thisFileName = MAIN_DIRECTORY+'/' + (str(timelist[t])).replace(' ', '_') + ceUniqueID +'.txt'
            # cloudElementsTextFile = open(thisFileName,'w')
            # -------------------------------------------------
            # ------ NETCDF File stuff for brightness temp stuff ------------------------------------
            thisFileName = userVariables.DIRS['mainDirStr'] + '/MERGnetcdfCEs/cloudElements' + \
                (str(timelist[t])).replace(' ', '_') + ceUniqueID + '.nc'
            currNetCDFCEData = Dataset(thisFileName, 'w', format='NETCDF4')
            currNetCDFCEData.description = 'Cloud Element ' + ceUniqueID + ' temperature data'
            currNetCDFCEData.calendar = 'standard'
            currNetCDFCEData.conventions = 'COARDS'
            # Dimensions
            currNetCDFCEData.createDimension('time', None)
            currNetCDFCEData.createDimension('lat', len(LAT[:, 0]))
            currNetCDFCEData.createDimension('lon', len(LON[0, :]))
            # Variables
            tempDims = ('time', 'lat', 'lon',)
            times = currNetCDFCEData.createVariable('time', 'f8', ('time',))
            times.units = 'hours since ' + str(timelist[t])[:-6]
            latitudes = currNetCDFCEData.createVariable('latitude', 'f8', ('lat',))
            longitudes = currNetCDFCEData.createVariable('longitude', 'f8', ('lon',))
            brightnesstemp = currNetCDFCEData.createVariable('brightnesstemp', 'f8', tempDims)
            brightnesstemp.units = 'Kelvin'
            # NETCDF data
            dates = [timelist[t] + timedelta(hours=0)]
            times[:] = date2num(dates, units=times.units)
            longitudes[:] = LON[0, :]
            longitudes.units = 'degrees_east'
            longitudes.long_name = 'Longitude'

            latitudes[:] = LAT[:, 0]
            latitudes.units = 'degrees_north'
            latitudes.long_name = 'Latitude'
            # -----------End most of NETCDF file stuff ------------------------------------

            # Generate array of zeros for brightness temperature
            brightnesstemp1 = ma.zeros((1, len(latitudes), len(longitudes))).astype('f8')
            finalCETRMMvalues = ma.zeros(brightnesstemp.shape)
            # Populate cloudElementLatLons by unpacking the original values from loc to get the actual value for lat and lon
            cloudElementNonZeros = cloudElement.nonzero()
            cloudyAreas = np.transpose(cloudElementNonZeros)
            cloudElementLatLons = np.zeros((cloudyAreas.shape[0], 3))
            cloudElementLatLons[:, 0] = cloudElementLat[cloudyAreas[:, 0]]
            cloudElementLatLons[:, 1] = cloudElementLon[cloudyAreas[:, 1]]
            cloudElementLatLons[:, 2] = cloudElement[cloudElementNonZeros]
            cloudElementLatLons = cloudElementLatLons.tolist()
            cloudElementLatLons = [tuple(l)for l in cloudElementLatLons]

            # This replaces the loop computation of brightnesstemp1, commented lines are a test
            brightnesstemp1[0, loc[0], loc[1]] = cloudElement
            brightnesstemp[:] = brightnesstemp1[:]
            currNetCDFCEData.close()

            # If other dataset (TRMM) assumed to be a precipitation dataset was entered
            if TRMMdirName:
                # ------ NETCDF File info for writing TRMM CE rainfall ------------------------------------
                thisFileName = userVariables.DIRS['mainDirStr']+'/TRMMnetcdfCEs/TRMM' + (str(timelist[t])).replace(' ', '_') + ceUniqueID + '.nc'
                currNetCDFTRMMData = Dataset(thisFileName, 'w', format='NETCDF4')
                currNetCDFTRMMData.description = 'Cloud Element ' + ceUniqueID + ' precipitation data'
                currNetCDFTRMMData.calendar = 'standard'
                currNetCDFTRMMData.conventions = 'COARDS'
                # Dimensions
                currNetCDFTRMMData.createDimension('time', None)
                currNetCDFTRMMData.createDimension('lat', len(LAT[:, 0]))
                currNetCDFTRMMData.createDimension('lon', len(LON[0, :]))

                # Variables
                TRMMprecip = ('time', 'lat', 'lon',)
                times = currNetCDFTRMMData.createVariable('time', 'f8', ('time',))
                times.units = 'hours since ' + str(timelist[t])[:-6]
                latitude = currNetCDFTRMMData.createVariable('latitude', 'f8', ('lat',))
                longitude = currNetCDFTRMMData.createVariable('longitude', 'f8', ('lon',))
                rainFallacc = currNetCDFTRMMData.createVariable('precipitation_Accumulation', 'f8', TRMMprecip)
                rainFallacc.units = 'mm'

                longitude[:] = LON[0, :]
                longitude.units = 'degrees_east'
                longitude.long_name = 'Longitude'

                latitude[:] = LAT[:, 0]
                latitude.units = 'degrees_north'
                latitude.long_name = 'Latitude'
                # -----------End most of NETCDF file stuff ------------------------------------
                chunkToInsert = np.copy(regriddedTRMM[loc[0], loc[1]])
                # chunkToInsert = regriddedTRMMCopy[loc[0],loc[1]]
                chunkToInsert[cloudElement == 0] = 0
                finalCETRMMvalues[0, loc[0], loc[1]] = chunkToInsert

                ceTRMMList = np.zeros((cloudyAreas.shape[0], 3))
                ceTRMMList[:, 0] = cloudElementLat[cloudyAreas[:, 0]]
                ceTRMMList[:, 1] = cloudElementLon[cloudyAreas[:, 1]]
                ceTRMMList[:, 2] = finalCETRMMvalues[0, np.floor(ceTRMMList[:, 0]).astype(int), np.floor(ceTRMMList[:, 1]).astype(int)]
                ceTRMMList = ceTRMMList.tolist()
                ceTRMMList = [tuple(l)for l in ceTRMMList]

                # Calculate the total precip associated with the feature
                precipTotal = np.sum(finalCETRMMvalues)
                rainFallacc[:] = finalCETRMMvalues[:]
                currNetCDFTRMMData.close()

                TRMMnumOfBoxes = np.count_nonzero(finalCETRMMvalues)
                TRMMArea = TRMMnumOfBoxes * userVariables.XRES * userVariables.YRES

                try:
                    maxCEprecipRate = np.max(finalCETRMMvalues[np.nonzero(finalCETRMMvalues)])
                    minCEprecipRate = np.min(finalCETRMMvalues[np.nonzero(finalCETRMMvalues)])
                except:
                    pass

            # Sort cloudElementLatLons by lats
            cloudElementLatLons.sort(key=lambda tup: tup[0])
            # Determine if the cloud element the shape
            cloudElementEpsilon = eccentricity(cloudElement)
            latCenter, lonCenter = ndimage.measurements.center_of_mass(cloudElement, labels=labels)
            # latCenter and lonCenter are given according to the particular array defining this CE
            # so you need to convert this value to the overall domain truth
            latCenter = cloudElementLat[round(latCenter)]
            lonCenter = cloudElementLon[round(lonCenter)]

            # Create the latLonBox
            latLonBox.append(min(cloudElementLon))
            latLonBox.append(min(cloudElementLat))
            latLonBox.append(max(cloudElementLon))
            latLonBox.append(max(cloudElementLat))

            nonEmptyLons = sum(sum(cloudElement) > 0)
            nonEmptyLats = sum(sum(cloudElement.transpose()) > 0)

            shape = max(nonEmptyLats, nonEmptyLons)

            cf = (((ndimage.minimum(cloudElement,
                    labels=labels)) / float((ndimage.maximum(cloudElement, labels=labels)))) * 100.0)
            tmin = ndimage.minimum(cloudElement, labels=labels)*1.
            tmax = ndimage.maximum(cloudElement, labels=labels)*1.
            avgBT = ndimage.mean(cloudElement, labels=labels)
            varBT = ndimage.variance(cloudElement, labels=labels)
            cloudElementCenter.append(latCenter)
            cloudElementCenter.append(lonCenter)
            # Populate the dictionary
            if TRMMdirName:
                cloudElementDict = {'uniqueID': ceUniqueID, 'cloudElementTime': timelist[t],
                    'cloudElementLatLon': cloudElementLatLons, 'cloudElementCenter': cloudElementCenter,
                    'cloudElementArea': cloudElementArea, 'cloudElementEccentricity': cloudElementEpsilon,
                    'cloudElementTmax': tmax, 'cloudElementTmin': tmin, 'cloudElementPrecipTotal': precipTotal,
                    'cloudElementLatLonTRMM': ceTRMMList, 'TRMMArea': TRMMArea, 'CETRMMmax': maxCEprecipRate,
                    'CETRMMmin': minCEprecipRate, 'cloudElementLatLonBox': latLonBox, 'cloudElementCF': cf,
                    'cloudElementBTavg': avgBT, 'cloudElementBTvar': varBT}
            else:
                cloudElementDict = {'uniqueID': ceUniqueID, 'cloudElementTime': timelist[t],
                    'cloudElementLatLon': cloudElementLatLons, 'cloudElementCenter': cloudElementCenter,
                    'cloudElementArea': cloudElementArea, 'cloudElementEccentricity': cloudElementEpsilon,
                    'cloudElementTmax': tmax, 'cloudElementTmin': tmin, 'cloudElementLatLonBox': latLonBox,
                    'cloudElementCF': cf, 'cloudElementBTavg': avgBT, 'cloudElementBTvar': varBT}

            # Data to be returned to parent function
            allCloudElementDicts.append(cloudElementDict)
        else:
            # Ensure only the non-zero elements are considered
            # Store intel in allCE file
            labels, _ = ndimage.label(cloudElement)
            cloudElementsFileString += ('\n-----------------------------------------------')
            cloudElementsFileString += ('\n\nTime is: %s' % (str(timelist[t])))
            cloudElementsFileString += ('\n\nceUniqueID is: %s' % ('F' + str(frameNum) + 'CE' + str(00)))

            cloudElementNonZeros = cloudElement.nonzero()
            cloudyAreas = np.transpose(cloudElementNonZeros)
            cloudElementLatLons = np.zeros((cloudyAreas.shape[0], 3))
            cloudElementLatLons[:, 0] = cloudElementLat[cloudyAreas[:, 0]]
            cloudElementLatLons[:, 1] = cloudElementLon[cloudyAreas[:, 1]]
            cloudElementLatLons[:, 2] = cloudElement[cloudElementNonZeros]
            cloudElementLatLons = cloudElementLatLons.tolist()
            cloudElementLatLons = [tuple(l)for l in cloudElementLatLons]

            cloudElementsFileString += ('\nLocation of rejected CE (lat,lon) points are: %s' % cloudElementLatLons)
            # latCenter and lonCenter are given according to the particular array defining this CE
            # so you need to convert this value to the overall domain truth
            latCenter, lonCenter = ndimage.measurements.center_of_mass(cloudElement, labels=labels)
            latCenter = cloudElementLat[round(latCenter)]
            lonCenter = cloudElementLon[round(lonCenter)]
            cloudElementsFileString += ('\nCenter (lat,lon) is: %.2f\t%.2f' % (latCenter, lonCenter))
            cloudElementsFileString += ('\nNumber of boxes are: %d' % numOfBoxes)
            cloudElementsFileString += ('\nArea is: %.4f km^2' % (cloudElementArea, ))
            cloudElementsFileString += ('\nAverage brightness temperature is: %.4f K' % ndimage.mean(cloudElement,
                labels=labels))
            cloudElementsFileString += ('\nMin brightness temperature is: %.4f K' % ndimage.minimum(cloudElement,
                labels=labels))
            cloudElementsFileString += ('\nMax brightness temperature is: %.4f K' % ndimage.maximum(cloudElement,
                labels=labels))
            cloudElementsFileString += ('\nBrightness temperature variance is: %.4f K'
                % ndimage.variance(cloudElement, labels=labels))
            cloudElementsFileString += ('\nConvective fraction is: %.4f ' % (((ndimage.minimum(cloudElement,
                labels=labels))/float((ndimage.maximum(cloudElement, labels=labels))))*100.0))
            cloudElementsFileString += ('\nEccentricity is: %.4f ' % (cloudElementEpsilon))
            cloudElementsFileString += ('\n-----------------------------------------------')

        # Reset list for the next CE
        cloudElementCenter = []
        cloudElement = []
        cloudElementLat = []
        cloudElementLon = []
        cloudElementLatLons = []
        brightnesstemp1 = []
        brightnesstemp = []
        finalCETRMMvalues = []
        ceTRMMList = []
        precipTotal = 0.0
        precip = []
        latLonBox = []

    return [allCloudElementDicts, cloudElementsFileString]
# **********************************************************************************************************************


def find_precip_rate(TRMMdirName, timelist, userVariables, graphVariables):
    counter = 0
    '''
    Purpose:: Determines the precipitation rates for MCSs found if TRMMdirName was not entered in
        find_cloud_elements this can be used

    Inputs:: TRMMdirName: a string representing the directory for the original TRMM netCDF files
        timelist: a list of python datatimes

    Returns:: allCEnodesTRMMdata: a list of dictionary of the TRMM data
        NB: also creates netCDF with TRMM data for each CE (for post processing) index
            in MAIN_DIRECTORY/TRMMnetcdfCEs

    Assumptions:: Assumes that find_cloud_elements was run without the TRMMdirName value

    '''
    allCEnodesTRMMdata = []
    precipTotal = 0.0

    os.chdir((userVariables.DIRS['mainDirStr'] + '/MERGnetcdfCEs/'))
    temporalRes = 3  # 3 hours for TRMM

    # Sort files
    files = filter(os.path.isfile, glob.glob('*.nc'))
    files.sort(key=lambda x: os.path.getmtime(x))

    for afile in files:
        fullFname = os.path.splitext(afile)[0]
        noFrameExtension = (fullFname.replace('_', '')).split('F')[0]
        ceUniqueID = 'F' + (fullFname.replace('_', '')).split('F')[1]
        fileDateTimeChar = (noFrameExtension.replace(':', '')).split('s')[1]
        fileDateTime = fileDateTimeChar.replace('-', '')
        fileDate = fileDateTime[:-6]
        fileHr1 = fileDateTime[-6:-4]

        cloudElementData = Dataset(afile, 'r', format='NETCDF4')
        brightnesstemp1 = cloudElementData.variables['brightnesstemp'][:, :, :]

        brightnesstemp = np.squeeze(brightnesstemp1, axis=0)

        if int(fileHr1) % temporalRes == 0:
            fileHr = fileHr1
        else:
            fileHr = (int(fileHr1) / temporalRes) * temporalRes

        if fileHr < 10:
            fileHr = '0' + str(fileHr)
        else:
            str(fileHr)

        TRMMfileName = TRMMdirName + '/3B42.' + str(fileDate) + '.'+str(fileHr) + '.7A.nc'
        TRMMData = Dataset(TRMMfileName, 'r', format='NETCDF4')
        precipRate = TRMMData.variables['pcp'][:, :, :]
        latsrawTRMMData = TRMMData.variables['latitude'][:]
        lonsrawTRMMData = TRMMData.variables['longitude'][:]
        lonsrawTRMMData[lonsrawTRMMData > 180] = lonsrawTRMMData[lonsrawTRMMData > 180] - 360.
        LONTRMM, LATTRMM = np.meshgrid(lonsrawTRMMData, latsrawTRMMData)

        # TODO: From this point onward 'LAT' and 'LON' are used but undefined.
        #       perhaps they have been replaced with LONTRMM and LATTRMM??
        nygrd = len(LAT[:, 0])
        nxgrd = len(LON[0, :])

        precipRateMasked = ma.masked_array(precipRate, mask=(precipRate < 0.0))
        # ---------regrid the TRMM data to the MERG dataset ----------------------------------
        # regrid using the do_regrid stuff from the Apache OCW
        regriddedTRMM = ma.zeros((0, nygrd, nxgrd))
        regriddedTRMM = utils.do_regrid(precipRateMasked[0, :, :], LATTRMM, LONTRMM, LAT, LON, order=1, mdi=-999999999)
        # ----------------------------------------------------------------------------------

        TRMMData.close()
        # ------ NETCDF File stuff ------------------------------------
        thisFileName = userVariables.DIRS['mainDirStr']+'/TRMMnetcdfCEs/' + fileDateTime + ceUniqueID + '.nc'
        currNetCDFTRMMData = Dataset(thisFileName, 'w', format='NETCDF4')
        currNetCDFTRMMData.description = 'Cloud Element ' + ceUniqueID + ' rainfall data'
        currNetCDFTRMMData.calendar = 'standard'
        currNetCDFTRMMData.conventions = 'COARDS'
        # dimensions
        currNetCDFTRMMData.createDimension('time', None)
        currNetCDFTRMMData.createDimension('lat', len(LAT[:, 0]))
        currNetCDFTRMMData.createDimension('lon', len(LON[0, :]))
        # variables
        TRMMprecip = ('time', 'lat', 'lon',)
        times = currNetCDFTRMMData.createVariable('time', 'f8', ('time',))
        times.units = 'hours since ' + fileDateTime[:-6]
        latitude = currNetCDFTRMMData.createVariable('latitude', 'f8', ('lat',))
        longitude = currNetCDFTRMMData.createVariable('longitude', 'f8', ('lon',))
        rainFallacc = currNetCDFTRMMData.createVariable('precipitation_Accumulation', 'f8', TRMMprecip)
        rainFallacc.units = 'mm'

        longitude[:] = LON[0, :]
        longitude.units = 'degrees_east'
        longitude.long_name = 'Longitude'

        latitude[:] = LAT[:, 0]
        latitude.units = 'degrees_north'
        latitude.long_name = 'Latitude'
        # -----------End most of NETCDF file stuff ------------------------------------

        finalCETRMMvalues = ma.zeros(brightnesstemp1.shape)
        # This block replaces the above commented loop
        finalCETRMMvalues[0, :, :] = regriddedTRMM
        finalCETRMMvalues[0, brightnesstemp <= 0] = 0

        precipTotal = np.sum(finalCETRMMvalues)
        rainFallacc[:] = finalCETRMMvalues
        currNetCDFTRMMData.close()
        TRMMnumOfBoxes = np.count_nonzero(finalCETRMMvalues)
        TRMMArea = TRMMnumOfBoxes * userVariables.XRES * userVariables.YRES

        try:
            minCEprecipRate = np.min(finalCETRMMvalues[np.nonzero(finalCETRMMvalues)])
        except:
            minCEprecipRate = 0.0

        try:
            maxCEprecipRate = np.max(finalCETRMMvalues[np.nonzero(finalCETRMMvalues)])
        except:
            maxCEprecipRate = 0.0

        # add info to CLOUDELEMENTSGRAPH
        # TODO try block
        for eachdict in graphVariables.CLOUD_ELEMENT_GRAPH.nodes(ceUniqueID):
            counter += 1
            if eachdict[1]['uniqueID'] == ceUniqueID:
                if 'cloudElementPrecipTotal' not in eachdict[1].keys():
                    eachdict[1]['cloudElementPrecipTotal'] = precipTotal
                if 'cloudElementLatLonTRMM' not in eachdict[1].keys():
                    eachdict[1]['cloudElementLatLonTRMM'] = finalCETRMMvalues
                if 'TRMMArea' not in eachdict[1].keys():
                    eachdict[1]['TRMMArea'] = TRMMArea
                if 'CETRMMmin' not in eachdict[1].keys():
                    eachdict[1]['CETRMMmin'] = minCEprecipRate
                if 'CETRMMmax' not in eachdict[1].keys():
                    eachdict[1]['CETRMMmax'] = maxCEprecipRate

        # Clean up
        precipTotal = 0.0
        latsrawTRMMData = []
        lonsrawTRMMData = []
        finalCETRMMvalues = []
        brightnesstemp = []

    return allCEnodesTRMMdata
# **********************************************************************************************************************


def find_cloud_clusters(CEGraph, userVariables, graphVariables):
    '''
    Purpose:: Determines the cloud clusters properties from the subgraphs in
        the graph i.e. prunes the graph according to the minimum depth

    Inputs:: CEGraph: a Networkx directed graph of the CEs with weighted edges
        according the area overlap between nodes (CEs) of consectuive frames

    Returns:: PRUNED_GRAPH: a Networkx directed graph of with CCs/ MCSs

    '''

    seenNode = []

    cloudClustersFile = open((userVariables.DIRS['mainDirStr'] + '/textFiles/cloudClusters.txt'), 'wb')

    for eachNode in CEGraph:
        # Check if the node has been seen before
        if eachNode not in dict(enumerate(zip(*seenNode))):
            # Look for all trees associated with node as the root
            thisPathDistanceAndLength = nx.single_source_dijkstra(CEGraph, eachNode)
            # Determine the actual shortestPath and minimum depth/length
            maxDepthAndMinPath = find_max_depth_and_min_path(thisPathDistanceAndLength, userVariables)
            if maxDepthAndMinPath:
                maxPathLength = maxDepthAndMinPath[0]
                shortestPath = maxDepthAndMinPath[1]

                # Add nodes and paths to PRUNED_GRAPH
                for i in xrange(len(shortestPath)):
                    if graphVariables.PRUNED_GRAPH.has_node(shortestPath[i]) is False:
                        graphVariables.PRUNED_GRAPH.add_node(shortestPath[i])

                    # Add edge if necessary
                    if i < (len(shortestPath)-1) and graphVariables.PRUNED_GRAPH.has_edge(shortestPath[i], shortestPath[i+1]) is False:
                        prunedGraphEdgeweight = CEGraph.get_edge_data(shortestPath[i], shortestPath[i+1])['weight']
                        graphVariables.PRUNED_GRAPH.add_edge(shortestPath[i], shortestPath[i+1], weight=prunedGraphEdgeweight)

                # Note information in a file for consideration later i.e. checking to see if it works
                cloudClustersFile.write('\nSubtree pathlength is %d and path is %s' % (maxPathLength, shortestPath))
                # update seenNode info
                seenNode.append(shortestPath)

    print 'pruned graph'
    print 'number of nodes are: ', graphVariables.PRUNED_GRAPH.number_of_nodes()
    print 'number of edges are: ', graphVariables.PRUNED_GRAPH.number_of_edges()
    print ('*'*80)

    # graphTitle = 'Cloud Clusters observed over somewhere during sometime'
    # plotting.draw_graph(graphVariables.PRUNED_GRAPH, graphTitle, userVariables.DIRS['mainDirStr'], graphVariables.edgeWeight)
    cloudClustersFile.close()

    return graphVariables.PRUNED_GRAPH
# **********************************************************************************************************************


def find_MCC(prunedGraph, userVariables, graphVariables):
    '''
    Purpose:: Determines if subtree is a MCC according to Laurent et al 1998 criteria

    Inputs:: prunedGraph: a Networkx Graph representing the CCs

    Returns:: finalMCCList: a list of list of tuples representing a MCC

    Assumptions:
        frames are ordered and are equally distributed in time e.g. hrly satellite images

    '''

    MCCList = []
    #MCSList = []
    definiteMCC = []
    definiteMCS = []
    eachList = []
    eachMCCList = []
    #maturing = False
    #decaying = False
    fNode = ''
    lNode = ''
    removeList = []
    imgCount = 0
    imgTitle = ''

    # maxShieldNode = ''
    orderedPath = []
    treeTraversalList = []
    definiteMCCFlag = False
    unDirGraph = nx.Graph()
    aSubGraph = nx.DiGraph()
    #definiteMCSFlag = False

    # Connected_components is not available for DiGraph, so generate graph as undirected
    unDirGraph = graphVariables.PRUNED_GRAPH.to_undirected()
    subGraph = nx.connected_component_subgraphs(unDirGraph)

    # For each path in the subgraphs determined
    for path in subGraph:
        # definite is a subTree provided the duration is longer than 3 hours
        if len(path.nodes()) > userVariables.MIN_MCS_DURATION:
            orderedPath = path.nodes()
            orderedPath.sort(key=lambda item: (len(item.split('C')[0]), item.split('C')[0]))
            # definiteMCS.append(orderedPath)

            # Build back DiGraph for checking purposes/paper purposes
            aSubGraph.add_nodes_from(path.nodes())
            for eachNode in path.nodes():
                if prunedGraph.predecessors(eachNode):
                    for node in prunedGraph.predecessors(eachNode):
                        aSubGraph.add_edge(node, eachNode, weight=graphVariables.edgeWeight[0])

                if prunedGraph.successors(eachNode):
                    for node in prunedGraph.successors(eachNode):
                        aSubGraph.add_edge(eachNode, node, weight=graphVariables.edgeWeight[0])
            imgTitle = 'CC'+str(imgCount+1)
            # Plotting.draw_graph(aSubGraph, imgTitle, MAIN_DIRECTORY, edgeWeight) #for eachNode in path:
            imgCount += 1
            # ----------end build back ---------------------------------------------
            mergeList, splitList = has_merges_or_splits(path, graphVariables)
            # Add node behavior regarding neutral, merge, split or both
            for node in path:
                if node in mergeList and node in splitList:
                    add_node_behavior_identifier(node, 'B', graphVariables)
                elif node in mergeList and node not in splitList:
                    add_node_behavior_identifier(node, 'M', graphVariables)
                elif node in splitList and node not in mergeList:
                    add_node_behavior_identifier(node, 'S', graphVariables)
                else:
                    add_node_behavior_identifier(node, 'N', graphVariables)
            # Do the first part of checking for the MCC feature
            # find the path
            treeTraversalList = traverse_tree(aSubGraph, orderedPath[0], [], [])
            # print 'treeTraversalList is ', treeTraversalList
            # check the nodes to determine if a MCC on just the area criteria (consecutive nodes meeting the area and temp requirements)
            MCCList = checked_nodes_MCC(prunedGraph, treeTraversalList, userVariables, graphVariables)
            for aDict in MCCList:
                for eachNode in aDict['fullMCSMCC']:
                    add_node_MCS_identifier(eachNode[0], eachNode[1], graphVariables)

            # o check for if MCCs overlap
            if MCCList:
                if len(MCCList) > 1:
                    for count in range(len(MCCList)):  # for eachDict in MCCList:
                        # if there are more than two lists
                        if count >= 1:
                            # and the first node in this list
                            eachList = list(x[0] for x in MCCList[count]['possMCCList'])
                            eachList.sort(key=lambda nodeID: (len(nodeID.split('C')[0]), nodeID.split('C')[0]))
                            if eachList:
                                fNode = eachList[0]
                                # Get the lastNode in the previous possMCC list
                                eachList = list(x[0] for x in MCCList[(count-1)]['possMCCList'])
                                eachList.sort(key=lambda nodeID: (len(nodeID.split('C')[0]), nodeID.split('C')[0]))
                                if eachList:
                                    lNode = eachList[-1]
                                    if lNode in graphVariables.CLOUD_ELEMENT_GRAPH.predecessors(fNode):
                                        for aNode in graphVariables.CLOUD_ELEMENT_GRAPH.predecessors(fNode):
                                            if aNode in eachList and aNode == lNode:
                                            # if edge_data is equal or less than to the exisitng edge in the tree append one to the other
                                                if graphVariables.CLOUD_ELEMENT_GRAPH.get_edge_data(aNode, fNode)['weight'] <= \
                                                        graphVariables.CLOUD_ELEMENT_GRAPH.get_edge_data(lNode, fNode)['weight']:
                                                    MCCList[count-1]['possMCCList'].extend(MCCList[count]['possMCCList'])
                                                    MCCList[count-1]['fullMCSMCC'].extend(MCCList[count]['fullMCSMCC'])
                                                    MCCList[count-1]['durationAandB'] += MCCList[count]['durationAandB']
                                                    MCCList[count-1]['CounterCriteriaA'] += MCCList[count]['CounterCriteriaA']
                                                    MCCList[count-1]['highestMCCnode'] = MCCList[count]['highestMCCnode']
                                                    MCCList[count-1]['frameNum'] = MCCList[count]['frameNum']
                                                    removeList.append(count)
                # Update the MCCList
                if removeList:
                    for i in removeList:
                        if (len(MCCList) - 1) > i:
                            del MCCList[i]
                            removeList = []

            # Check if the nodes also meet the duration criteria and the shape crieria
            for eachDict in MCCList:
                # Order the fullMCSMCC list, then run maximum extent and eccentricity criteria
                if (eachDict['durationAandB'] * userVariables.TRES) >= userVariables.MINIMUM_DURATION and \
                        (eachDict['durationAandB'] * userVariables.TRES) <= userVariables.MAXIMUM_DURATION:
                    eachList = list(x[0] for x in eachDict['fullMCSMCC'])
                    eachList.sort(key=lambda nodeID: (len(nodeID.split('C')[0]), nodeID.split('C')[0]))
                    eachMCCList = list(x[0] for x in eachDict['possMCCList'])
                    eachMCCList.sort(key=lambda nodeID: (len(nodeID.split('C')[0]), nodeID.split('C')[0]))

                    # Update the nodemcsidentifer behavior
                    # Find the first element eachMCCList in eachList, and ensure everything ahead of it is indicated as 'I',
                    # Find last element in eachMCCList in eachList and ensure everything after it is indicated as 'D'
                    # Ensure that everything between is listed as 'M'
                    for eachNode in eachList[:(eachList.index(eachMCCList[0]))]:
                        add_node_MCS_identifier(eachNode, 'I', graphVariables)

                    add_node_MCS_identifier(eachMCCList[0], 'M', graphVariables)

                    for eachNode in eachList[(eachList.index(eachMCCList[-1])+1):]:
                        add_node_MCS_identifier(eachNode, 'D', graphVariables)

                    # update definiteMCS list
                    for eachNode in orderedPath[(orderedPath.index(eachMCCList[-1])+1):]:
                        add_node_MCS_identifier(eachNode, 'D', graphVariables)

                    # run maximum extent and eccentricity criteria
                    _, definiteMCCFlag = max_extent_and_eccentricity(eachList, userVariables, graphVariables)
                    # maxExtentNode, definiteMCCFlag = max_extent_and_eccentricity(eachList)
                    # print 'maxExtentNode, definiteMCCFlag ', maxExtentNode, definiteMCCFlag
                    if definiteMCCFlag is True:
                        definiteMCC.append(eachList)

            definiteMCS.append(orderedPath)

            # Reset for next subGraph
            aSubGraph.clear()
            orderedPath = []
            MCCList = []
            #MCSList =[]
            #definiteMCSFlag = False

    return definiteMCC, definiteMCS
# **********************************************************************************************************************


def traverse_tree(subGraph, node, stack, checkedNodes=None):
    '''
    Purpose:: To traverse a tree using a modified depth-first iterative deepening (DFID) search algorithm

    Input:: subGraph: a Networkx DiGraph representing a CC
            lengthOfsubGraph: an integer representing the length of the subgraph
            node: a string representing the node currently being checked
            stack: a list of strings representing a list of nodes in a stack functionality
                    i.e. Last-In-First-Out (LIFO) for sorting the information from each visited node
            checkedNodes: a list of strings representing the list of the nodes in the traversal

    Returns:: checkedNodes: a list of strings representing the list of the nodes in the traversal

    Assumptions: frames are ordered and are equally distributed in time e.g. hrly satellite images

    '''

    if len(checkedNodes) == len(subGraph):
        return checkedNodes

    if not checkedNodes:
        stack = []
        checkedNodes.append(node)

    # Check one level infront first...if something does exisit, stick it at the front of the stack
    upOneLevel = subGraph.predecessors(node)
    downOneLevel = subGraph.successors(node)
    for parent in upOneLevel:
        if parent not in checkedNodes and parent not in stack:
            for child in downOneLevel:
                if child not in checkedNodes and child not in stack:
                    stack.insert(0, child)

            stack.insert(0, parent)

    for child in downOneLevel:
        if child not in checkedNodes and child not in stack:
            if len(subGraph.predecessors(child)) > 1 or node in checkedNodes:
                stack.insert(0, child)
            else:
                stack.append(child)

    for eachNode in stack:
        if eachNode not in checkedNodes:
            checkedNodes.append(eachNode)
            return traverse_tree(subGraph, eachNode, stack, checkedNodes)

    return checkedNodes
# **********************************************************************************************************************


def checked_nodes_MCC(prunedGraph, nodeList, userVariables, graphVariables):
    '''
    Purpose ::Determine if this path is (or is part of) a MCC and provides
        preliminary information regarding the stages of the feature

    Inputs:: prunedGraph: a Networkx Graph representing all the cloud clusters
        nodeList: list of strings (CE ID) from the traversal

    Returns:: potentialMCCList: list of dictionaries representing all possible MCC within the path
            dictionary = {'possMCCList':[(node,'I')], 'fullMCSMCC':[(node,'I')], 'CounterCriteriaA': CounterCriteriaA,
                         'durationAandB': durationAandB}
    '''

    counterCriteriaAFlag = False
    counterCriteriaBFlag = False
    initiationFlag = False
    maturityFlag = False
    decayFlag = False
    thisdict = {}  # Will have the same items as the cloudElementDict
    cloudElementAreaB = 0.0
    # cloudElementAreaA = 0.0
    # epsilon = 0.0
    # frameNum =0
    oldNode = ''
    potentialMCCList = []
    # durationAandB = 0

    # Check for if the list contains only one string/node
    if type(nodeList) is str:
        oldNode = nodeList
        nodeList = []
        nodeList.append(oldNode)

    for node in nodeList:
        thisdict = this_dict(node, graphVariables)
        counterCriteriaAFlag = False
        counterCriteriaBFlag = False
        # existingFrameFlag = False

        if thisdict['cloudElementArea'] >= userVariables.OUTER_CLOUD_SHIELD_AREA:
            counterCriteriaAFlag = True
            initiationFlag = True
            maturityFlag = False

            # Check if criteriaA is met
            cloudElementAreaA, _ = check_criteria(thisdict['cloudElementLatLon'], userVariables.OUTER_CLOUD_SHIELD_TEMPERATURE, userVariables)
            # cloudElementAreaA, criteriaA = check_criteria(thisdict['cloudElementLatLon'], OUTER_CLOUD_SHIELD_TEMPERATURE)
            # TODO: calculate the eccentricity at this point and read over????or create a new field in the dict

            if cloudElementAreaA >= userVariables.OUTER_CLOUD_SHIELD_AREA:
                # Check if criteriaB is met
                cloudElementAreaB, criteriaB = check_criteria(thisdict['cloudElementLatLon'], userVariables.INNER_CLOUD_SHIELD_TEMPERATURE, userVariables)

                # If Criteria A and B have been met, then the MCC is initiated, i.e. store node as potentialMCC
                if cloudElementAreaB >= userVariables.INNER_CLOUD_SHIELD_AREA:
                    # TODO: add another field to the dictionary for the OUTER_AREA_SHIELD area
                    counterCriteriaBFlag = True
                    # Append this information on to the dictionary
                    add_info_this_dict(node, cloudElementAreaB, criteriaB, graphVariables)
                    initiationFlag = False
                    maturityFlag = True
                    stage = 'M'
                    potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,
                                                       counterCriteriaAFlag, counterCriteriaBFlag)
                else:
                    # Criteria B failed
                    counterCriteriaBFlag = False
                    if initiationFlag is True:
                        stage = 'I'
                        potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,
                                                           counterCriteriaAFlag, counterCriteriaBFlag)

                    elif (initiationFlag is False and maturityFlag is True) or decayFlag is True:
                        decayFlag = True
                        maturityFlag = False
                        stage = 'D'
                        potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,
                                                           counterCriteriaAFlag, counterCriteriaBFlag)
            else:
                # Criteria A failed
                counterCriteriaAFlag = False
                counterCriteriaBFlag = False
                # Add as a CE before or after the main feature
                if initiationFlag is True or (initiationFlag is False and maturityFlag is True):
                    stage = 'I'
                    potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,
                                                       counterCriteriaAFlag, counterCriteriaBFlag)
                elif (initiationFlag is False and maturityFlag is False) or decayFlag is True:
                    stage = 'D'
                    decayFlag = True
                    potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,
                                                       counterCriteriaAFlag, counterCriteriaBFlag)
                elif (initiationFlag is False and maturityFlag is False and decayFlag is False):
                    stage = 'I'
                    potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,
                                                       counterCriteriaAFlag, counterCriteriaBFlag)
        else:
            # Criteria A failed
            counterCriteriaAFlag = False
            counterCriteriaBFlag = False
            # Add as a CE before or after the main feature
            if initiationFlag is True or (initiationFlag is False and maturityFlag is True):
                stage = 'I'
                potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,
                                                   counterCriteriaAFlag, counterCriteriaBFlag)
            elif (initiationFlag is False and maturityFlag is False) or decayFlag is True:
                stage = 'D'
                decayFlag = True
                potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,
                                                   counterCriteriaAFlag, counterCriteriaBFlag)
            elif (initiationFlag is False and maturityFlag is False and decayFlag is False):
                stage = 'I'
                potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,
                                                   counterCriteriaAFlag, counterCriteriaBFlag)
    return potentialMCCList
# **********************************************************************************************************************


def update_MCC_list(prunedGraph, potentialMCCList, node, stage, counterCriteriaAFlag, counterCriteriaBFlag):
    '''
    Purpose::
        Utility function to determine if a path is (or is part of) a MCC and provides
               preliminary information regarding the stages of the feature

    Input::
        prunedGraph: a Networkx Graph representing all the cloud clusters
        potentialMCCList: a list of dictionaries representing the possible MCCs within a path
        node: a string representing the cloud element currently being assessed
        counterCriteriaAFlag: a boolean value indicating whether the node meets the MCC criteria A according to Laurent et al
        counterCriteriaBFlag: a boolean value indicating whether the node meets the MCC criteria B according to Laurent et al

    Output::
        potentialMCCList: list of dictionaries representing all possible MCC within the path
             dictionary = {'possMCCList':[(node,'I')], 'fullMCSMCC':[(node,'I')], 'CounterCriteriaA': CounterCriteriaA, 'durationAandB': durationAandB}

    '''

    existingFrameFlag = False
    #existingMCSFrameFlag = False
    predecessorsFlag = False
    predecessorsMCSFlag = False
    successorsFlag = False
    successorsMCSFlag = False
    frameNum = 0

    frameNum = int((node.split('CE')[0]).split('F')[1])
    if potentialMCCList == []:
        # List empty
        stage = 'I'
        if counterCriteriaAFlag is True and counterCriteriaBFlag is True:
            potentialMCCList.append({'possMCCList': [(node, stage)], 'fullMCSMCC': [(node, stage)], 'CounterCriteriaA': 1,
                 'durationAandB': 1, 'highestMCCnode': node, 'frameNum': frameNum})
        elif counterCriteriaAFlag is True and counterCriteriaBFlag is False:
            potentialMCCList.append({'possMCCList': [], 'fullMCSMCC': [(node, stage)], 'CounterCriteriaA': 1,
                 'durationAandB': 0, 'highestMCCnode': '', 'frameNum': 0})
        elif counterCriteriaAFlag is False and counterCriteriaBFlag is False:
            potentialMCCList.append({'possMCCList': [], 'fullMCSMCC': [(node, stage)], 'CounterCriteriaA': 0,
                 'durationAandB': 0, 'highestMCCnode': '', 'frameNum': 0})
    else:
        # List not empty
        predecessorsFlag, index = is_there_a_link(prunedGraph, 1, node, potentialMCCList, 1)
        if predecessorsFlag is True:
            for eachNode in potentialMCCList[index]['possMCCList']:
                if int((eachNode[0].split('CE')[0]).split('F')[1]) == frameNum:
                    existingFrameFlag = True
            # This MUST come after the check for the existing frame
            if counterCriteriaAFlag is True and counterCriteriaBFlag is True:
                stage = 'M'
                potentialMCCList[index]['possMCCList'].append((node, stage))
                potentialMCCList[index]['fullMCSMCC'].append((node, stage))

            if existingFrameFlag is False:
                if counterCriteriaAFlag is True and counterCriteriaBFlag is True:
                    stage = 'M'
                    potentialMCCList[index]['CounterCriteriaA'] += 1
                    potentialMCCList[index]['durationAandB'] += 1
                    if frameNum > potentialMCCList[index]['frameNum']:
                        potentialMCCList[index]['frameNum'] = frameNum
                        potentialMCCList[index]['highestMCCnode'] = node
                    return potentialMCCList

                # If this frameNum doesn't exist and this frameNum is less than the MCC node max frame Num (including 0), then append to fullMCSMCC list
                if frameNum > potentialMCCList[index]['frameNum'] or potentialMCCList[index]['frameNum'] == 0:
                    stage = 'I'
                    if counterCriteriaAFlag is True and counterCriteriaBFlag is False:
                        potentialMCCList.append({'possMCCList': [], 'fullMCSMCC': [(node, stage)], 'CounterCriteriaA': 1,
                             'durationAandB': 0, 'highestMCCnode': '', 'frameNum': 0})
                        return potentialMCCList
                    elif counterCriteriaAFlag is False and counterCriteriaBFlag is False:
                        potentialMCCList.append({'possMCCList': [], 'fullMCSMCC': [(node, stage)], 'CounterCriteriaA': 0,
                             'durationAandB': 0, 'highestMCCnode': '', 'frameNum': 0})
                        return potentialMCCList

            # If predecessor and this frame number already exist in the MCC list, add the current node to the fullMCSMCC list
            if existingFrameFlag is True:
                if counterCriteriaAFlag is True and counterCriteriaBFlag is False:
                    potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                    potentialMCCList[index]['CounterCriteriaA'] += 1
                    return potentialMCCList
                if counterCriteriaAFlag is False:
                    potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                    return potentialMCCList

        if predecessorsFlag is False:
            successorsFlag, index = is_there_a_link(prunedGraph, 2, node, potentialMCCList, 2)

            if successorsFlag is True:
                for eachNode in potentialMCCList[index]['possMCCList']:
                    if int((eachNode[0].split('CE')[0]).split('F')[1]) == frameNum:
                        existingFrameFlag = True

                if counterCriteriaAFlag is True and counterCriteriaBFlag is True:
                    stage = 'M'
                    potentialMCCList[index]['possMCCList'].append((node, stage))
                    potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                    if frameNum > potentialMCCList[index]['frameNum'] or potentialMCCList[index]['frameNum'] == 0:
                        potentialMCCList[index]['frameNum'] = frameNum
                        potentialMCCList[index]['highestMCCnode'] = node
                    return potentialMCCList

                if existingFrameFlag is False:
                    if stage == 'M':
                        stage = 'D'
                    if counterCriteriaAFlag is True and counterCriteriaBFlag is True:
                        potentialMCCList[index]['CounterCriteriaA'] += 1
                        potentialMCCList[index]['durationAandB'] += 1
                    elif counterCriteriaAFlag is True:
                        potentialMCCList[index]['CounterCriteriaA'] += 1
                    elif counterCriteriaAFlag is False:
                        potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                        return potentialMCCList
    # If predecessor and this frame number already exist in the MCC list, add the current node to the fullMCSMCC list
                else:
                    if counterCriteriaAFlag is True and counterCriteriaBFlag is False:
                        potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                        potentialMCCList[index]['CounterCriteriaA'] += 1
                        return potentialMCCList
                    if counterCriteriaAFlag is False:
                        potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                        return potentialMCCList

        # If this node isn't connected to existing MCCs check if it is connected to existing MCSs ...
        if predecessorsFlag is False and successorsFlag is False:
            stage = 'I'
            predecessorsMCSFlag, index = is_there_a_link(prunedGraph, 1, node, potentialMCCList, 2)
            if predecessorsMCSFlag is True:
                if counterCriteriaAFlag is True and counterCriteriaBFlag is True:
                    potentialMCCList[index]['possMCCList'].append((node, 'M'))
                    potentialMCCList[index]['fullMCSMCC'].append((node, 'M'))
                    potentialMCCList[index]['durationAandB'] += 1
                    if frameNum > potentialMCCList[index]['frameNum']:
                        potentialMCCList[index]['frameNum'] = frameNum
                        potentialMCCList[index]['highestMCCnode'] = node
                    return potentialMCCList

                if potentialMCCList[index]['frameNum'] == 0 or frameNum <= potentialMCCList[index]['frameNum']:
                    if counterCriteriaAFlag is True and counterCriteriaBFlag is False:
                        potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                        potentialMCCList[index]['CounterCriteriaA'] += 1
                        return potentialMCCList
                    elif counterCriteriaAFlag is False:
                        potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                        return potentialMCCList
            else:
                successorsMCSFlag, index = is_there_a_link(prunedGraph, 2, node, potentialMCCList, 2)
                if successorsMCSFlag is True:
                    if counterCriteriaAFlag is True and counterCriteriaBFlag is True:
                        potentialMCCList[index]['possMCCList'].append((node, 'M'))
                        potentialMCCList[index]['fullMCSMCC'].append((node, 'M'))
                        potentialMCCList[index]['durationAandB'] += 1
                        if frameNum > potentialMCCList[index]['frameNum']:
                            potentialMCCList[index]['frameNum'] = frameNum
                            potentialMCCList[index]['highestMCCnode'] = node
                        return potentialMCCList

                    if potentialMCCList[index]['frameNum'] == 0 or frameNum <= potentialMCCList[index]['frameNum']:
                        if counterCriteriaAFlag is True and counterCriteriaBFlag is False:
                            potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                            potentialMCCList[index]['CounterCriteriaA'] += 1
                            return potentialMCCList
                        elif counterCriteriaAFlag is False:
                            potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                            return potentialMCCList

            # If this node isn't connected to existing MCCs or MCSs, create a new one ...
            if predecessorsFlag is False and predecessorsMCSFlag is False and successorsFlag is False and\
               successorsMCSFlag is False:
                if counterCriteriaAFlag is True and counterCriteriaBFlag is True:
                    potentialMCCList.append({'possMCCList': [(node, stage)], 'fullMCSMCC': [(node, stage)],
                         'CounterCriteriaA': 1, 'durationAandB': 1, 'highestMCCnode': node, 'frameNum': frameNum})
                elif counterCriteriaAFlag is True and counterCriteriaBFlag is False:
                    potentialMCCList.append({'possMCCList': [], 'fullMCSMCC': [(node, stage)], 'CounterCriteriaA': 1,
                         'durationAandB': 0, 'highestMCCnode': '', 'frameNum': 0})
                elif counterCriteriaAFlag is False and counterCriteriaBFlag is False:
                    potentialMCCList.append({'possMCCList': [], 'fullMCSMCC': [(node, stage)], 'CounterCriteriaA': 0,
                         'durationAandB': 0, 'highestMCCnode': '', 'frameNum': 0})

    return potentialMCCList
# **********************************************************************************************************************


def is_there_a_link(prunedGraph, upOrDown, node, potentialMCCList, whichList):
    '''
    Purpose::
        Utility script for update_MCC_list mostly because there is no Pythonic way to break out of nested loops

    Input::
        prunedGraph:a Networkx Graph representing all the cloud clusters
        upOrDown: an integer representing 1- to do predecesor check and 2 - to do successor checked_nodes_MCC
        node: a string representing the cloud element currently being assessed
        potentialMCCList: a list of dictionaries representing the possible MCCs within a path
        whichList: an integer representing which list ot check in the dictionary; 1- possMCCList, 2- fullMCSMCC

    Output::
        thisFlag: a boolean representing whether the list passed has in the parent or child of the node
        index: an integer representing the location in the potentialMCCList where thisFlag occurs

    '''

    thisFlag = False
    index = -1
    checkList = ''
    if whichList == 1:
        checkList = 'possMCCList'
    elif whichList == 2:
        checkList = 'fullMCSMCC'

    # Check parents
    if upOrDown == 1:
        for aNode in prunedGraph.predecessors(node):
            # Reset the index counter for this node search through potentialMCCList
            index = -1
            for mccDict in potentialMCCList:
                index += 1
                if aNode in list(x[0] for x in mccDict[checkList]):
                    thisFlag = True
                    # Get out of looping so as to avoid the flag being written over when another node in the predecessor
                    # list is checked
                    return thisFlag, index

    # Check children
    if upOrDown == 2:
        for aNode in prunedGraph.successors(node):
            # Reset the index counter for this node search through potentialMCCList
            index = -1
            for mccDict in potentialMCCList:
                index += 1

                if aNode in list(x[0] for x in mccDict[checkList]):
                    thisFlag = True
                    return thisFlag, index

    return thisFlag, index
# **********************************************************************************************************************


def max_extent_and_eccentricity(eachList, userVariables, graphVariables):
    '''
    Purpose::
        Perform the final check for MCC based on maximum extent and eccentricity criteria

    Input::
        eachList: a list of strings  representing the node of the possible MCCs within a path

    Output::
        maxShieldNode: a string representing the node with the maximum maxShieldNode
        definiteMCCFlag: a boolean indicating that the MCC has met all requirements

    '''

    maxShieldNode = ''
    maxShieldArea = 0.0
    # maxShieldEccentricity = 0.0
    definiteMCCFlag = False

    if eachList:
        for eachNode in eachList:
            if (this_dict(eachNode, graphVariables)['nodeMCSIdentifier'] == 'M' or this_dict(eachNode, graphVariables)['nodeMCSIdentifier'] == 'D') and\
                    this_dict(eachNode, graphVariables)['cloudElementArea'] > maxShieldArea:
                maxShieldNode = eachNode
                maxShieldArea = this_dict(eachNode, graphVariables)['cloudElementArea']

        # maxShieldEccentricity = this_dict(maxShieldNode)['cloudElementEccentricity']
        if this_dict(maxShieldNode, graphVariables)['cloudElementEccentricity'] >= userVariables.ECCENTRICITY_THRESHOLD_MIN and \
                this_dict(maxShieldNode, graphVariables)['cloudElementEccentricity'] <= userVariables.ECCENTRICITY_THRESHOLD_MAX:
            # Criteria met
            definiteMCCFlag = True

    return maxShieldNode, definiteMCCFlag
#  **********************************************************************************************************************


def find_max_depth_and_min_path(thisPathDistanceAndLength, userVariables):
    '''
    Purpose::
        To determine the maximum depth and min path for the headnode

    Input::
        tuple of dictionaries representing the shortest distance and paths for a node in the tree as 
        returned by nx.single_source_dijkstra
        thisPathDistanceAndLength({distance}, {path})
            {distance} = nodeAsString, valueAsInt, {path} = nodeAsString, pathAsList

    Output::
        tuple of the max pathLength and min pathDistance as a tuple (like what was input)
            minDistanceAndMaxPath = ({distance},{path})
    '''

    maxPathLength = 0
    minPath = 0

    # maxPathLength for the node in question
    maxPathLength = max(len(values) for values in thisPathDistanceAndLength[1].values())

    # If the duration is shorter then the min MCS length, then don't store!
    if maxPathLength < userVariables.MIN_MCS_DURATION:  # MINIMUM_DURATION:
        minDistanceAndMaxPath = ()

    # else find the min path and max depth
    else:
        # Max path distance for the node in question
        minPath = max(values for values in thisPathDistanceAndLength[0].values())

        # Check to determine the shortest path from the longest paths returned
        for pathDistance, path in itertools.izip(thisPathDistanceAndLength[0].values(), thisPathDistanceAndLength[1].values()):
            pathLength = len(path)
            # If pathLength is the same as the maxPathLength, then look the pathDistance to determine if the min
            if pathLength == maxPathLength:
                if pathDistance <= minPath:
                    minPath = pathLength
                    # Store details if absolute minPath and deepest
                    minDistanceAndMaxPath = (pathDistance, path)
    return minDistanceAndMaxPath
# **********************************************************************************************************************


def this_dict(thisNode, graphVariables):
    '''
    Purpose::
        Return dictionary from graph if node exist in tree

    Input::
        thisNode: a string representing the CE to get the information for

    Output ::
        eachdict[1]: a dictionary representing the info associated with thisNode from the graph

    '''

    for eachdict in graphVariables.CLOUD_ELEMENT_GRAPH.nodes(thisNode):
        if eachdict[1]['uniqueID'] == thisNode:
            return eachdict[1]
# **********************************************************************************************************************


def check_criteria(thisCloudElementLatLon, aTemperature, userVariables):
    '''
    Purpose:: Determine if criteria B is met for a CEGraph

    Inputs:: thisCloudElementLatLon: an array of (lat,lon, t_bb)
        aTemperature:a integer representing the temperature maximum for masking

    Returns:: cloudElementArea: a floating-point number representing the area in the array that meet the criteria 

    '''

    cloudElementCriteriaBLatLon = []

    _, ceCounter = ndimage.measurements.label(thisCloudElementLatLon, structure=userVariables.STRUCTURING_ELEMENT)

    allCriteriaB = []

    # Determine min and max values in lat and lon, then use this to generate teh array from LAT,LON meshgrid

    minLat = min(x[0] for x in thisCloudElementLatLon)
    maxLat = max(x[0] for x in thisCloudElementLatLon)
    minLon = min(x[1] for x in thisCloudElementLatLon)
    maxLon = max(x[1] for x in thisCloudElementLatLon)

    minLatIndex = np.argmax(LAT[:, 0] == minLat)
    maxLatIndex = np.argmax(LAT[:, 0] == maxLat)
    minLonIndex = np.argmax(LON[0, :] == minLon)
    maxLonIndex = np.argmax(LON[0, :] == maxLon)

    criteriaBframe = ma.zeros(((abs(maxLatIndex - minLatIndex)+1), (abs(maxLonIndex - minLonIndex)+1)))

    for x in thisCloudElementLatLon:
        # To store the values of the subset in the new array, remove the minLatIndex and minLonindex from the
        # index given in the original array to get the indices for the new array
        criteriaBframe[(np.argmax(LAT[:, 0] == x[0]) - minLatIndex), (np.argmax(LON[0, :] == x[1]) - minLonIndex)] = x[2]

    # Keep only those values < aTemperature
    tempMask = ma.masked_array(criteriaBframe, mask=(criteriaBframe >= aTemperature), fill_value=0)

    # Get the actual values that the mask returned
    criteriaB = ma.zeros(criteriaBframe.shape).astype('f8')
    criteriaB[~tempMask.mask] = tempMask[~tempMask.mask]

    for _ in xrange(ceCounter):
        # [0] is time dimension. Determine the actual values from the data
        # loc is a masked array
        # ***** returns elements down then across thus (6,4) is 6 arrays deep of size 4
        try:
            loc = ndimage.find_objects(criteriaB)[0]
        except:
            # This would mean that no objects were found meeting criteria B
            cloudElementArea = 0.0
            continue

        try:
            cloudElementCriteriaB = ma.zeros(criteriaB.shape)
            cloudElementCriteriaB = criteriaB[loc]
        except:
            print 'YIKESS'
            print 'ceCounter ', ceCounter, criteriaB.shape
            print 'criteriaB ', criteriaB

        cloudElementCriteriaBNonZeros = cloudElementCriteriaB.nonzero()
        cloudElement = np.transpose(cloudElementCriteriaBNonZeros)
        cloudElementCriteriaBLatLon = np.zeros((cloudElement.shape[0], 3))
        cloudElementCriteriaBLatLon[:, 0] = LAT[cloudElement[:, 1], 0]
        cloudElementCriteriaBLatLon[:, 1] = LON[0, cloudElement[:, 2]]
        cloudElementCriteriaBLatLon[:, 2] = cloudElementCriteriaB[cloudElementCriteriaBNonZeros]
        cloudElementCriteriaBLatLon = cloudElementCriteriaBLatLon.tolist()
        cloudElementCriteriaBLatLon = [tuple(l)for l in cloudElementCriteriaBLatLon]

        cloudElementArea = np.count_nonzero(cloudElementCriteriaB) * userVariables.XRES * userVariables.YRES

        tempMask = []
        criteriaB = []
        cloudElementCriteriaB = []

        allCriteriaB.append((cloudElementArea, cloudElementCriteriaBLatLon))

    return max(allCriteriaB, key=lambda x: x[0]) if allCriteriaB != [] else cloudElementArea, cloudElementCriteriaBLatLon
# **********************************************************************************************************************


def has_merges_or_splits(nodeList, graphVariables):
    '''
    Purpose:: Determine if nodes within a path defined from shortest_path splittingNodeDict
    Inputs:: nodeList: list of strings representing the nodes from a path
    Returns:: splitList: a list of strings representing all the nodes in the path that split
        mergeList: a list of strings representing all the nodes in the path that merged
    '''

    mergeList = []
    splitList = []

    for node, numParents in graphVariables.PRUNED_GRAPH.in_degree(nodeList).items():
        if numParents > 1:
            mergeList.append(node)

    for node, numChildren in graphVariables.PRUNED_GRAPH.out_degree(nodeList).items():
        if numChildren > 1:
            splitList.append(node)
    # Sort
    splitList.sort(key=lambda item: (len(item.split('C')[0]), item.split('C')[0]))
    mergeList.sort(key=lambda item: (len(item.split('C')[0]), item.split('C')[0]))

    return mergeList, splitList
# **********************************************************************************************************************


def all_ancestors(path, aNode):
    '''
    Purpose:: Utility script to provide the path leading up to a nodeList

    Inputs:: path: a list of strings representing the nodes in the path
        aNode: a string representing a node to be checked for parents

    Returns:: path: a list of strings representing the list of the nodes connected to aNode through its parents
        numOfChildren: an integer representing the number of parents of the node passed
    '''

    numOfParents = PRUNED_GRAPH.in_degree(aNode)
    try:
        if PRUNED_GRAPH.predecessors(aNode) and numOfParents <= 1:
            path = path + PRUNED_GRAPH.predecessors(aNode)
            thisNode = PRUNED_GRAPH.predecessors(aNode)[0]
            return all_ancestors(path, thisNode)
        else:
            path = path + aNode
            return path, numOfParents
    except:
        return path, numOfParents
# **********************************************************************************************************************


def all_descendants(path, aNode):
    '''
    Purpose:: Utility script to provide the path leading up to a nodeList

    Inputs:: path: a list of strings representing the nodes in the path
        aNode: a string representing a node to be checked for children

    Returns:: path: a list of strings representing the list of the nodes connected to aNode through its children
        numOfChildren: an integer representing the number of children of the node passed
    '''

    numOfChildren = PRUNED_GRAPH.out_degree(aNode)
    try:
        if PRUNED_GRAPH.successors(aNode) and numOfChildren <= 1:
            path = path + PRUNED_GRAPH.successors(aNode)
            thisNode = PRUNED_GRAPH.successors(aNode)[0]
            return all_descendants(path, thisNode)
        else:
            path = path + aNode
            # i.e. PRUNED_GRAPH.predecessors(aNode) is empty
            return path, numOfChildren
    except:
        # i.e. PRUNED_GRAPH.predecessors(aNode) threw an exception
        return path, numOfChildren
# **********************************************************************************************************************


def add_info_this_dict(thisNode, cloudElementArea, criteriaB, graphVariables):
    '''
    Purpose:: Update original dictionary node with information

    Inputs:: thisNode: a string representing the unique ID of a node
        cloudElementArea: a floating-point number representing the area of the cloud element
        criteriaB: a masked array of floating-point numbers representing the lat,lons meeting the criteria

    Returns:: None

    '''

    for eachdict in graphVariables.CLOUD_ELEMENT_GRAPH.nodes(thisNode):
        if eachdict[1]['uniqueID'] == thisNode:
            eachdict[1]['CriteriaBArea'] = cloudElementArea
            eachdict[1]['CriteriaBLatLon'] = criteriaB
    return
# **********************************************************************************************************************


def add_node_behavior_identifier(thisNode, nodeBehaviorIdentifier, graphVariables):
    '''
    Purpose:: add an identifier to the node dictionary to indicate splitting, merging or neither node

    Inputs:: thisNode: a string representing the unique ID of a node
        nodeBehaviorIdentifier: a string representing the behavior S- split, M- merge, B- both split and merge, 
        N- neither split or merge

    Returns:: None

    '''
    for eachdict in graphVariables.CLOUD_ELEMENT_GRAPH.nodes(thisNode):
        if eachdict[1]['uniqueID'] == thisNode:
            if 'nodeBehaviorIdentifier' not in eachdict[1].keys():
                eachdict[1]['nodeBehaviorIdentifier'] = nodeBehaviorIdentifier
    return
# **********************************************************************************************************************


def add_node_MCS_identifier(thisNode, nodeMCSIdentifier, graphVariables):
    '''
    Purpose:: Add an identifier to the node dictionary to indicate splitting, merging or neither node

    Inputs:: thisNode: a string representing the unique ID of a node
        nodeMCSIdentifier: a string representing the stage of the MCS lifecycle  'I' for Initiation, 'M' for Maturity,
        'D' for Decay

    Returns:: None

    '''

    for eachdict in graphVariables.CLOUD_ELEMENT_GRAPH.nodes(thisNode):
        if eachdict[1]['uniqueID'] == thisNode:
            if 'nodeMCSIdentifier' not in eachdict[1].keys():
                eachdict[1]['nodeMCSIdentifier'] = nodeMCSIdentifier
    return
# **********************************************************************************************************************


def update_node_MCS_identifier(thisNode, nodeMCSIdentifier, graphVariables):
    '''
    Purpose:: Update an identifier to the node dictionary to indicate splitting, merging or neither node

    Inputs:: thisNode: thisNode: a string representing the unique ID of a node
        nodeMCSIdentifier: a string representing the stage of the MCS lifecyle  'I' for Initiation, 'M' for Maturity, 'D' for Decay

    Returns:: None

    '''
    for eachdict in graphVariables.CLOUD_ELEMENT_GRAPH.nodes(thisNode):
        if eachdict[1]['uniqueID'] == thisNode:
            eachdict[1]['nodeMCSIdentifier'] = nodeMCSIdentifier

    return
# **********************************************************************************************************************


def eccentricity(cloudElementLatLon):
    '''
    Purpose:: Determines the eccentricity (shape) of contiguous boxes
        Values tending to 1 are more circular by definition, whereas
        values tending to 0 are more linear

    Inputs:: cloudElementLatLon: 2D array in (lat,lon) representing T_bb contiguous squares

    Returns:: epsilon: a floating-point representing the eccentricity of the matrix passed

    '''

    epsilon = 0.0

    # nonEmptyLons = sum(sum(cloudElementLatLon) > 0)
    # nonEmptyLats = sum(sum(cloudElementLatLon.transpose()) > 0)

    # I think this is what's wanted
    sh = cloudElementLatLon.shape
    nonEmptyLons = sh[1]
    nonEmptyLats = sh[0]

    lonEigenvalues = 1.0 * nonEmptyLats / (nonEmptyLons + 0.001)  # for long oval on y axis
    latEigenvalues = 1.0 * nonEmptyLons / (nonEmptyLats + 0.001)  # for long oval on x-axs
    epsilon = min(latEigenvalues, lonEigenvalues)

    #     THIS LOOP APPEARS TO BE UNNECESSARY
    # Loop over all lons and determine longest (non-zero) col
    # Loop over all lats and determine longest (non-zero) row
    # sh = cloudElementLatLon.shape

#    for _ in cloudElementLatLon:
#        assign a matrix to determine the legit values
#
#        nonEmptyLons = sum(sum(cloudElementLatLon) > 0)
#        nonEmptyLats = sum(sum(cloudElementLatLon.transpose()) > 0)
#        assert_sameval(nonEmptyLons,sh[1],nonEmptyLats,sh[0])

#        lonEigenvalues = 1.0 * nonEmptyLats / (nonEmptyLons + 0.001) #for long oval on y axis
#        latEigenvalues = 1.0 * nonEmptyLons / (nonEmptyLats + 0.001) #for long oval on x-axs
#        epsilon = min(latEigenvalues, lonEigenvalues)

    # Assert(epsilon_0==epsilon)
    return epsilon


def assert_sameval(n1, n2, n3, n4):
    if(n1 != n2 or n3 != n4):
        bad = 1
    return
# **********************************************************************************************************************


def cloud_element_overlap(currentCELatLons, previousCELatLons, xRes, yRes):
    '''
    Purpose::
        Determines the percentage overlap between two list of lat-lons passed

    Input::
        currentCELatLons: a list of tuples for the current CE
        previousCELatLons: a list of tuples for the other CE being considered
        xRes: a floating-point represening the distance between longitudes (in km)
        yRes: a floating-point representing the distance between latitudes (in km)

    Output::
        percentageOverlap: a floating-point representing the number of overlapping lat_lon tuples
        areaOverlap: a floating-point number representing the area overlapping

    '''

    latlonprev = []
    latloncurr = []
    count = 0
    percentageOverlap = 0.0
    areaOverlap = 0.0

    # Remove the temperature from the tuples for currentCELatLons and previousCELatLons then check for overlap
    latlonprev = [(x[0], x[1]) for x in previousCELatLons]
    latloncurr = [(x[0], x[1]) for x in currentCELatLons]

    # Find overlap
    count = len(list(set(latloncurr) & set(latlonprev)))

    # Find area overlap
    areaOverlap = count * xRes * yRes

    # Find percentage
    percentageOverlap = max(((count * 1.0) / (len(latloncurr) * 1.0)), ((count * 1.0) / (len(latlonprev) * 1.0)))

    return percentageOverlap, areaOverlap
# **********************************************************************************************************************
