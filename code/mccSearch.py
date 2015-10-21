from datetime import timedelta
import glob
import itertools
from netCDF4 import Dataset, date2num
import numpy as np
import numpy.ma as ma
import os
from scipy import ndimage
import json

import networkx as nx

#GTG modules
import utils
import plotting


#----------------------- GLOBAL VARIABLES --------------------------
# --------------------- User defined variables ---------------------
#FYI the lat lon values are not necessarily inclusive of the points given. These are the limits
#the first point closest the the value (for the min) from the MERG data is used, etc.
LATMIN = '5.0' #min latitude; -ve values in the SH e.g. 5S = -5
LATMAX = '19.0' #max latitude; -ve values in the SH e.g. 5S = -5 20.0
LONMIN = '-5.0' #min longitude; -ve values in the WH e.g. 59.8W = -59.8 -30
LONMAX = '9.0' #min longitude; -ve values in the WH e.g. 59.8W = -59.8  30
XRES = 4.0              #x direction spatial resolution in km
YRES = 4.0              #y direction spatial resolution in km
TRES = 1                #temporal resolution in hrs
LAT_DISTANCE = 111.0    #the avg distance in km for 1deg lat for the region being considered
LON_DISTANCE = 111.0    #the avg distance in km for 1deg lon for the region being considered
STRUCTURING_ELEMENT = [[0, 1, 0],
                       [1, 1, 1],
                       [0, 1, 0]
                      ] #the matrix for determining the pattern for the contiguous boxes and must
                        #have same rank of the matrix it is being compared against
                        #criteria for determining cloud elements and edges
T_BB_MAX = 243  #warmest temp to allow (-30C to -55C according to Morel and Sensi 2002)
T_BB_MIN = 218  #cooler temp for the center of the system
CONVECTIVE_FRACTION = 0.90 #the min temp/max temp that would be expected in a CE.. this is highly
                           #conservative (only a 10K difference)
MIN_MCS_DURATION = 3    #minimum time for a MCS to exist
AREA_MIN = 2400.0       #minimum area for CE criteria in km^2 according to Vila et al. (2008) is 2400
MIN_OVERLAP = 10000.00   #km^2  from Williams and Houze 1987, indir ref in Arnaud et al 1992

#---the MCC criteria
ECCENTRICITY_THRESHOLD_MAX = 1.0  #tending to 1 is a circle e.g. hurricane,
ECCENTRICITY_THRESHOLD_MIN = 0.70 #tending to 0 is a linear e.g. squall line
OUTER_CLOUD_SHIELD_AREA = 80000.0 #km^2
INNER_CLOUD_SHIELD_AREA = 30000.0 #km^2
OUTER_CLOUD_SHIELD_TEMPERATURE = 233 #in K
INNER_CLOUD_SHIELD_TEMPERATURE = 213 #in K
MINIMUM_DURATION = 6  #min number of frames the MCC must exist for (assuming hrly frames, MCCs is 6hrs)
MAXIMUM_DURATION = 24 #max number of framce the MCC can last for
#------------------- End user defined Variables -------------------
edgeWeight = [1, 2, 3] #weights for the graph edges
#graph object fo the CEs meeting the criteria
CLOUD_ELEMENT_GRAPH = nx.DiGraph()
#graph meeting the CC criteria
PRUNED_GRAPH = nx.DiGraph()
#get lat lons from iomethods.py

#------------------------ End GLOBAL VARS -------------------------
#************************ Begin Functions *************************
#**********************************************************************************************************************
def find_cloud_elements(mergImgs, timelist, mainStrDir, lat, lon, TRMMdirName=None):
    '''
    Purpose:: Determines the contiguous boxes for a given time of the satellite images i.e. each frame
        using scipy ndimage package

    Inputs:: mergImgs: masked numpy array in (time,lat,lon),T_bb representing the satellite data. This is masked based 
        on the maximum acceptable temperature, T_BB_MAX
        timelist: a list of python datatimes
        mainStrDir: a string representing the path to main directory where data will be written
        lat: a 2D array of the set of latitudes from the files opened
        lon: a 2D array of the set of longitudes from the files opened
        TRMMdirName (optional): string representing the path where to find the TRMM datafiles

    Returns::
        CLOUD_ELEMENT_GRAPH: a Networkx directed graph where each node contains the information in cloudElementDict
        The nodes are determined according to the area of contiguous squares. The nodes are linked through weighted
        edges.

        cloudElementDict = {'uniqueID': unique tag for this CE,
                            'cloudElementTime': time of the CE,
                            'cloudElementLatLon': (lat,lon,value) of MERG data of CE,
                            'cloudElementCenter':list of floating-point [lat,lon] representing the CE's center
                            'cloudElementArea':floating-point representing the area of the CE,
                            'cloudElementEccentricity': floating-point representing the shape of the CE,
                            'cloudElementTmax':integer representing the maximum Tb in CE,
                            'cloudElementTmin': integer representing the minimum Tb in CE,
                            'cloudElementPrecipTotal':floating-point representing the sum of all rainfall in CE if
                                                      TRMMdirName entered,
                            'cloudElementLatLonTRMM':(lat,lon,value) of TRMM data in CE if TRMMdirName entered,
                            'TRMMArea': floating-point representing the CE if TRMMdirName entered,
                            'CETRMMmax':floating-point representing the max rate in the CE if TRMMdirName entered,
                            'CETRMMmin':floating-point representing the min rate in the CE if TRMMdirName entered}
    Assumptions::
        Assumes we are dealing with MERG data which is 4kmx4km resolved, thus the smallest value
        required according to Vila et al. (2008) is 2400km^2
        therefore, 2400/16 = 150 contiguous squares
    '''

    global MAIN_DIRECTORY
    MAIN_DIRECTORY = mainStrDir
    global LAT
    LAT = lat
    global LON
    LON = lon

    cloudElementsJSON = []  #list of the key, value objects associated with a CE in the graph
    edges = []     #list of the nodes connected to a given CE
    latLonBox = [] #list of extreme points for the min fitting box around a CE [lon_min, lat_min, lon_max, lat_max]
    shape = 0      #max(num of non-zero boxes in lat, num of non-zero boxes in lon) for the CE
    cf = 0.0       #convective fraction 
    tmin = 0.0     #IR Tmin for the CE
    tmax = 0.0     #IR Tmax for the CE

    frame = ma.empty((1, mergImgs.shape[1], mergImgs.shape[2]))
    ceCounter = 0
    frameceCounter = 0
    frameNum = 0
    cloudElementEpsilon = 0.0
    cloudElementDict = {}
    cloudElementCenter = []     #list with two elements [lat,lon] for the center of a CE
    prevFrameCEs = []           #list for CEs in previous frame
    currFrameCEs = []           #list for CEs in current frame
    cloudElementLat = []        #list for a particular CE's lat values
    cloudElementLon = []        #list for a particular CE's lon values
    cloudElementLatLons = []    #list for a particular CE's (lat,lon) values

    TIR_min = 0.0
    TIR_max = 0.0
    temporalRes = 3 # TRMM data is 3 hourly
    precipTotal = 0.0
    ceTRMMList = []
    precip = []

    nygrd = len(LAT[:, 0])
    nxgrd = len(LON[0, :])

    #openfile for storing ALL cloudElement information
    cloudElementsFile = open((MAIN_DIRECTORY + '/textFiles/cloudElements.txt'), 'wb')
    #openfile for storing cloudElement information meeting user criteria i.e. MCCs in this case
    cloudElementsUserFile = open((MAIN_DIRECTORY + '/textFiles/cloudElementsUserFile.txt'), 'w')


    filenameJSON = MAIN_DIRECTORY + '/textFiles/graphJSON.txt'
    #NB in the TRMM files the info is hours since the time thus 00Z file has in 01, 02 and 03 times
    for t in xrange(mergImgs.shape[0]):
        #-------------------------------------------------
        # #textfile name for saving the data for arcgis
        # thisFileName = MAIN_DIRECTORY+'/' + (str(timelist[t])).replace(' ', '_') + '.txt'
        # cloudElementsTextFile = open(thisFileName,'w')
        #-------------------------------------------------

        #determine contiguous locations with temeperature below the warmest temp i.e. cloudElements in each frame
        frame, ceCounter = ndimage.measurements.label(mergImgs[t, :, :], structure=STRUCTURING_ELEMENT)
        frameceCounter = 0
        frameNum += 1

        #for each of the areas identified, check to determine if it a valid CE via an area and T requirement
        for count in xrange(ceCounter):
            #[0] is time dimension. Determine the actual values from the data
            #loc is a masked array
            try:
                loc = ndimage.find_objects(frame == (count + 1))[0]
            except Exception, e:
                print 'Error is ', e
                continue


            cloudElement = mergImgs[t, :, :][loc]
            labels, _ = ndimage.label(cloudElement)

            #determine the true lats and lons for this particular CE
            cloudElementLat = LAT[loc[0], 0]
            cloudElementLon = LON[0, loc[1]]

            #determine number of boxes in this cloudelement
            numOfBoxes = np.count_nonzero(cloudElement)
            cloudElementArea = numOfBoxes * XRES * YRES

            #If the area is greater than the area required, or if the area is smaller than the suggested area,
            #check if it meets a convective fraction requirement consider as CE

            if cloudElementArea >= AREA_MIN or (cloudElementArea < AREA_MIN and \
                    ((ndimage.minimum(cloudElement, labels=labels)) / \
                        float((ndimage.maximum(cloudElement, labels=labels)))) < CONVECTIVE_FRACTION):
                #get some time information and labeling info
                frameceCounter += 1
                ceUniqueID = 'F' + str(frameNum) + 'CE' + str(frameceCounter)
                #-------------------------------------------------
                #textfile name for accesing CE data using MATLAB code
                # thisFileName = MAIN_DIRECTORY+'/' + (str(timelist[t])).replace(' ', '_') + ceUniqueID +'.txt'
                # cloudElementsTextFile = open(thisFileName,'w')
                #-------------------------------------------------
                # ------ NETCDF File stuff for brightness temp stuff ------------------------------------
                thisFileName = MAIN_DIRECTORY + '/MERGnetcdfCEs/cloudElements' + \
                                (str(timelist[t])).replace(' ', '_') + ceUniqueID + '.nc'
                currNetCDFCEData = Dataset(thisFileName, 'w', format='NETCDF4')
                currNetCDFCEData.description = 'Cloud Element '+ ceUniqueID + ' temperature data'
                currNetCDFCEData.calendar = 'standard'
                currNetCDFCEData.conventions = 'COARDS'
                # dimensions
                currNetCDFCEData.createDimension('time', None)
                currNetCDFCEData.createDimension('lat', len(LAT[:, 0]))
                currNetCDFCEData.createDimension('lon', len(LON[0, :]))
                # variables
                tempDims = ('time', 'lat', 'lon',)
                times = currNetCDFCEData.createVariable('time', 'f8', ('time',))
                times.units = 'hours since '+ str(timelist[t])[:-6]
                latitudes = currNetCDFCEData.createVariable('latitude', 'f8', ('lat',))
                longitudes = currNetCDFCEData.createVariable('longitude', 'f8', ('lon',))
                brightnesstemp = currNetCDFCEData.createVariable('brightnesstemp', 'i16', tempDims)
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

                #generate array of zeros for brightness temperature
                brightnesstemp1 = ma.zeros((1, len(latitudes), len(longitudes))).astype('int16')
                #-----------End most of NETCDF file stuff ------------------------------------

                #if other dataset (TRMM) assumed to be a precipitation dataset was entered
                if TRMMdirName:
                    #------------------TRMM stuff -------------------------------------------------
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

                    #open TRMM file for the resolution info and to create the appropriate sized grid
                    TRMMfileName = TRMMdirName + '/3B42.' + fileDate + '.' + str(fileHr) + '.7A.nc'

                    TRMMData = Dataset(TRMMfileName, 'r', format='NETCDF4')
                    precipRate = TRMMData.variables['pcp'][:, :, :]
                    latsrawTRMMData = TRMMData.variables['latitude'][:]
                    lonsrawTRMMData = TRMMData.variables['longitude'][:]
                    lonsrawTRMMData[lonsrawTRMMData > 180] = lonsrawTRMMData[lonsrawTRMMData > 180] - 360.
                    LONTRMM, LATTRMM = np.meshgrid(lonsrawTRMMData, latsrawTRMMData)

                    # nygrdTRMM = len(LATTRMM[:, 0])
                    # nxgrdTRMM = len(LONTRMM[0, :])
                    precipRateMasked = ma.masked_array(precipRate, mask=(precipRate < 0.0))
                    #---------regrid the TRMM data to the MERG dataset ----------------------------------
                    #regrid using the do_regrid stuff from the Apache OCW
                    regriddedTRMM = ma.zeros((0, nygrd, nxgrd))
                    regriddedTRMM = utils.do_regrid(precipRateMasked[0, :, :], LATTRMM, LONTRMM, LAT, LON, order=1,\
                                    mdi=-999999999)
                    #----------------------------------------------------------------------------------

                    # #get the lat/lon info from cloudElement
                    #get the lat/lon info from the file
                    #latCEStart = LAT[0][0]
                    #latCEEnd = LAT[-1][0]
                    #lonCEStart = LON[0][0]
                    #lonCEEnd = LON[0][-1]

                    #get the lat/lon info for TRMM data (different resolution)
                    # latStartT = utils.find_nearest(latsrawTRMMData, latCEStart)
                    # latEndT = utils.find_nearest(latsrawTRMMData, latCEEnd)
                    # lonStartT = utils.find_nearest(lonsrawTRMMData, lonCEStart)
                    # lonEndT = utils.find_nearest(lonsrawTRMMData, lonCEEnd)
                    # Unused since CEPrecipRate isn't used and these are just inputs
                    # latStartIndex = np.where(latsrawTRMMData == latStartT)
                    # latEndIndex = np.where(latsrawTRMMData == latEndT)
                    # lonStartIndex = np.where(lonsrawTRMMData == lonStartT)
                    # lonEndIndex = np.where(lonsrawTRMMData == lonEndT)

                    #get the relevant TRMM info
                    # Unused Variable
                    # CEprecipRate = precipRate[:, (latStartIndex[0][0]-1):latEndIndex[0][0], \
                    #        (lonStartIndex[0][0]-1):lonEndIndex[0][0]]
                    TRMMData.close()

                    # ------ NETCDF File info for writing TRMM CE rainfall ------------------------------------
                    thisFileName = MAIN_DIRECTORY+'/TRMMnetcdfCEs/TRMM' + (str(timelist[t])).replace(' ', '_') + ceUniqueID +'.nc'
                    currNetCDFTRMMData = Dataset(thisFileName, 'w', format='NETCDF4')
                    currNetCDFTRMMData.description = 'Cloud Element ' + ceUniqueID + ' precipitation data'
                    currNetCDFTRMMData.calendar = 'standard'
                    currNetCDFTRMMData.conventions = 'COARDS'
                    # dimensions
                    currNetCDFTRMMData.createDimension('time', None)
                    currNetCDFTRMMData.createDimension('lat', len(LAT[:,0]))
                    currNetCDFTRMMData.createDimension('lon', len(LON[0,:]))

                    # variables
                    TRMMprecip = ('time', 'lat', 'lon',)
                    times = currNetCDFTRMMData.createVariable('time', 'f8', ('time',))
                    times.units = 'hours since '+ str(timelist[t])[:-6]
                    latitude = currNetCDFTRMMData.createVariable('latitude', 'f8', ('lat',))
                    longitude = currNetCDFTRMMData.createVariable('longitude', 'f8', ('lon',))
                    rainFallacc = currNetCDFTRMMData.createVariable('precipitation_Accumulation', 'f8', TRMMprecip)
                    rainFallacc.units = 'mm'

                    longitude[:] = LON[0,:]
                    longitude.units = 'degrees_east'
                    longitude.long_name = 'Longitude'

                    latitude[:] = LAT[:,0]
                    latitude.units = 'degrees_north'
                    latitude.long_name = 'Latitude'

                    finalCETRMMvalues = ma.zeros((brightnesstemp.shape))
                    #-----------End most of NETCDF file stuff ------------------------------------

                #populate cloudElementLatLons by unpacking the original values from loc to get the actual value for lat and lon
                #TODO: KDW - too dirty... play with itertools.izip or zip and the enumerate with this
                #           as cloudElement is masked
                for index, value in np.ndenumerate(cloudElement):
                    if value != 0:
                        latIndex, lonIndex = index
                        latLonTuple = (cloudElementLat[latIndex], cloudElementLon[lonIndex], value)

                        #generate the comma separated file for GIS
                        cloudElementLatLons.append(latLonTuple)

                        #temp data for CE NETCDF file
                        brightnesstemp1[0, int(np.where(LAT[:,0] == cloudElementLat[latIndex])[0]), \
                                int(np.where(LON[0,:] == cloudElementLon[lonIndex])[0])] = value

                        if TRMMdirName:
                            finalCETRMMvalues[0, int(np.where(LAT[:,0] == cloudElementLat[latIndex])[0]), \
                                int(np.where(LON[0,:] == cloudElementLon[lonIndex])[0])] = \
                                regriddedTRMM[int(np.where(LAT[:,0] == cloudElementLat[latIndex])[0]), \
                                int(np.where(LON[0,:] == cloudElementLon[lonIndex])[0])]
                            ceTRMMList.append((cloudElementLat[latIndex], cloudElementLon[lonIndex], \
                                finalCETRMMvalues[0,cloudElementLat[latIndex], cloudElementLon[lonIndex]]))

                brightnesstemp[:] = brightnesstemp1[:]
                currNetCDFCEData.close()

                if TRMMdirName:
                    #calculate the total precip associated with the feature
                    for index, value in np.ndenumerate(finalCETRMMvalues):
                        precipTotal += value
                        precip.append(value)

                    rainFallacc[:] = finalCETRMMvalues[:]
                    currNetCDFTRMMData.close()
                    TRMMnumOfBoxes = np.count_nonzero(finalCETRMMvalues)
                    TRMMArea = TRMMnumOfBoxes * XRES * YRES
                    try:
                        maxCEprecipRate = np.max(finalCETRMMvalues[np.nonzero(finalCETRMMvalues)])
                        minCEprecipRate = np.min(finalCETRMMvalues[np.nonzero(finalCETRMMvalues)])
                    except:
                        pass

                #sort cloudElementLatLons by lats
                cloudElementLatLons.sort(key=lambda tup: tup[0])

                #determine if the cloud element the shape
                cloudElementEpsilon = eccentricity(cloudElement)
                cloudElementsUserFile.write('\n\nTime is: %s' %(str(timelist[t])))
                cloudElementsUserFile.write('\nceUniqueID is: %s' %ceUniqueID)
                latCenter, lonCenter = ndimage.measurements.center_of_mass(cloudElement, labels=labels)

                #latCenter and lonCenter are given according to the particular array defining this CE
                #so you need to convert this value to the overall domain truth
                latCenter = cloudElementLat[round(latCenter)]
                lonCenter = cloudElementLon[round(lonCenter)]

                #create the latLonBox
                latLonBox.append(min(cloudElementLon))
                latLonBox.append(min(cloudElementLat))
                latLonBox.append(max(cloudElementLon))
                latLonBox.append(max(cloudElementLat))

                for _ in cloudElement:
                    #assign a matrix to determine the legit values
                    nonEmptyLons = sum(sum(cloudElement) > 0)
                    nonEmptyLats = sum(sum(cloudElement.transpose()) > 0)

                shape = max(nonEmptyLats, nonEmptyLons)

                cloudElementsUserFile.write('\nCenter (lat,lon) is: %.2f\t%.2f' %(latCenter, lonCenter))
                cloudElementCenter.append(latCenter)
                cloudElementCenter.append(lonCenter)
                cloudElementsUserFile.write('\nNumber of boxes are: %d' %numOfBoxes)
                cloudElementsUserFile.write('\nArea is: %.4f km^2' %(cloudElementArea))
                cloudElementsUserFile.write('\nAverage brightness temperature is: %.4f K' %ndimage.mean(cloudElement, \
                    labels=labels))
                cloudElementsUserFile.write('\nMin brightness temperature is: %.4f K' %ndimage.minimum(cloudElement, \
                    labels=labels))
                cloudElementsUserFile.write('\nMax brightness temperature is: %.4f K' %ndimage.maximum(cloudElement, \
                    labels=labels))
                cloudElementsUserFile.write('\nBrightness temperature variance is: %.4f K' \
                    %ndimage.variance(cloudElement, labels=labels))
                cloudElementsUserFile.write('\nConvective fraction is: %.4f ' %(((ndimage.minimum(cloudElement, \
                    labels=labels)) / float((ndimage.maximum(cloudElement, labels=labels)))) * 100.0))
                cloudElementsUserFile.write('\nEccentricity is: %.4f ' %(cloudElementEpsilon))
                #populate the dictionary
                if TRMMdirName:
                    cloudElementDict = {'uniqueID': ceUniqueID, 'cloudElementTime': timelist[t],\
                        'cloudElementLatLon': cloudElementLatLons, 'cloudElementCenter':cloudElementCenter, \
                        'cloudElementArea':cloudElementArea, 'cloudElementEccentricity':cloudElementEpsilon, \
                        'cloudElementTmax':TIR_max, 'cloudElementTmin': TIR_min, 'cloudElementPrecipTotal':precipTotal,\
                        'cloudElementLatLonTRMM':ceTRMMList, 'TRMMArea': TRMMArea, 'CETRMMmax':maxCEprecipRate, \
                        'CETRMMmin':minCEprecipRate}

                else:
                    cloudElementDict = {'uniqueID': ceUniqueID, 'cloudElementTime': timelist[t],\
                        'cloudElementLatLon': cloudElementLatLons, 'cloudElementCenter':cloudElementCenter, \
                        'cloudElementArea':cloudElementArea, 'cloudElementEccentricity':cloudElementEpsilon, \
                        'cloudElementTmax':TIR_max, 'cloudElementTmin': TIR_min,}

                #current frame list of CEs
                currFrameCEs.append(cloudElementDict)

                #draw the graph node
                CLOUD_ELEMENT_GRAPH.add_node(ceUniqueID, cloudElementDict)

                if frameNum != 1:
                    for cloudElementDict in prevFrameCEs:
                        percentageOverlap, areaOverlap = cloud_element_overlap(cloudElementLatLons, cloudElementDict['cloudElementLatLon'])

                        #change weights to int as the built in shortest path chokes on floating pts according to Networkx doc
                        #according to Goyens et al, two CEs are considered related if there is atleast 95% overlap between 
                        #them for consecutive imgs a max of 2 hrs apart
                        if percentageOverlap >= 0.95:
                            CLOUD_ELEMENT_GRAPH.add_edge(cloudElementDict['uniqueID'], ceUniqueID, weight=edgeWeight[0])
                            edges.append(cloudElementDict['uniqueID'])

                        elif percentageOverlap >= 0.90 and percentageOverlap < 0.95 :
                            CLOUD_ELEMENT_GRAPH.add_edge(cloudElementDict['uniqueID'], ceUniqueID, weight=edgeWeight[1])
                            edges.append(cloudElementDict['uniqueID'])

                        elif areaOverlap >= MIN_OVERLAP:
                            CLOUD_ELEMENT_GRAPH.add_edge(cloudElementDict['uniqueID'], ceUniqueID, weight=edgeWeight[2])
                            edges.append(cloudElementDict['uniqueID'])

                # get some data for the JSON object which will only store the graph CEs and connected edges 
                cf = (((ndimage.minimum(cloudElement, \
                    labels=labels)) / float((ndimage.maximum(cloudElement, labels=labels)))) * 100.0)
                tmin = ndimage.minimum(cloudElement, labels=labels)*1.
                tmax = ndimage.maximum(cloudElement, labels=labels)*1.
                if edges:
                    cloudElementsJSON.append({'cloudElement': ceUniqueID, 'time': str(timelist[t]),\
                        'area':cloudElementArea, 'Tmax': tmax, 'Tmin': tmin,'center':cloudElementCenter,\
                        'convective_fraction': cf, 'lat_lon_box': latLonBox, 'shape': shape, 'edges':edges, 'eccentricity':cloudElementEpsilon})

            else:
                #TODO: remove this else as we only wish for the CE details
                #ensure only the non-zero elements are considered
                #store intel in allCE file
                labels, _ = ndimage.label(cloudElement)
                cloudElementsFile.write('\n-----------------------------------------------')
                cloudElementsFile.write('\n\nTime is: %s' %(str(timelist[t])))
                cloudElementsFile.write('\n\nceUniqueID is: %s' %('F' + str(frameNum) + 'CE' + str(00)))
                # cloudElementLat = LAT[loc[0],0]
                # cloudElementLon = LON[0,loc[1]]

                #populate cloudElementLatLons by unpacking the original values from loc
                #TODO: KDW - too dirty... play with itertools.izip or zip and the enumerate with this
                #           as cloudElement is masked
                for index, value in np.ndenumerate(cloudElement):
                    if value != 0 :
                        latIndex, lonIndex = index
                        latLonTuple = (cloudElementLat[latIndex], cloudElementLon[lonIndex])
                        cloudElementLatLons.append(latLonTuple)

                cloudElementsFile.write('\nLocation of rejected CE (lat,lon) points are: %s' % cloudElementLatLons)
                #latCenter and lonCenter are given according to the particular array defining this CE
                #so you need to convert this value to the overall domain truth
                latCenter, lonCenter = ndimage.measurements.center_of_mass(cloudElement, labels=labels)
                latCenter = cloudElementLat[round(latCenter)]
                lonCenter = cloudElementLon[round(lonCenter)]
                cloudElementsFile.write('\nCenter (lat,lon) is: %.2f\t%.2f' % (latCenter, lonCenter))
                cloudElementsFile.write('\nNumber of boxes are: %d' % numOfBoxes)
                cloudElementsFile.write('\nArea is: %.4f km^2' % (cloudElementArea, ))
                cloudElementsFile.write('\nAverage brightness temperature is: %.4f K' % ndimage.mean(cloudElement, \
                    labels=labels))
                cloudElementsFile.write('\nMin brightness temperature is: %.4f K' % ndimage.minimum(cloudElement, \
                    labels=labels))
                cloudElementsFile.write('\nMax brightness temperature is: %.4f K' % ndimage.maximum(cloudElement, \
                    labels=labels))
                cloudElementsFile.write('\nBrightness temperature variance is: %.4f K' \
                    % ndimage.variance(cloudElement, labels=labels))
                cloudElementsFile.write('\nConvective fraction is: %.4f ' % (((ndimage.minimum(cloudElement, \
                    labels=labels))/float((ndimage.maximum(cloudElement, labels=labels))))*100.0))
                cloudElementsFile.write('\nEccentricity is: %.4f ' % (cloudElementEpsilon))
                cloudElementsFile.write('\n-----------------------------------------------')

            #reset list for the next CE
            cloudElementCenter = []
            cloudElement = []
            cloudElementLat = []
            cloudElementLon = []
            cloudElementLatLons = []
            brightnesstemp1 = []
            brightnesstemp = []
            finalCETRMMvalues = []
            #CEprecipRate =[]
            ceTRMMList = []
            precipTotal = 0.0
            precip = []
            edges = []
            latLonBox = []


        #reset for the next time
        prevFrameCEs = []
        prevFrameCEs = currFrameCEs
        currFrameCEs = []


    cloudElementsFile.close()
    cloudElementsUserFile.close()
    #if using ARCGIS data store code, uncomment this file close line
    #cloudElementsTextFile.close

    #write JSON file
    with open(filenameJSON, 'w+') as f:
        json.dump(cloudElementsJSON,f)

    #clean up graph - remove parent and childless nodes
    outAndInDeg = CLOUD_ELEMENT_GRAPH.degree_iter()
    toRemove = [node[0] for node in outAndInDeg if node[1] < 1]
    CLOUD_ELEMENT_GRAPH.remove_nodes_from(toRemove)

    print 'number of nodes are: ', CLOUD_ELEMENT_GRAPH.number_of_nodes()
    print 'number of edges are: ', CLOUD_ELEMENT_GRAPH.number_of_edges()
    print ('*'*80)

    #hierachial graph output
    # graphTitle = 'Cloud Elements observed over somewhere from 0000Z to 0000Z'
    # plotting.draw_graph(CLOUD_ELEMENT_GRAPH, graphTitle, MAIN_DIRECTORY, edgeWeight)

    return CLOUD_ELEMENT_GRAPH
#**********************************************************************************************************************
def find_precip_rate(TRMMdirName, timelist):
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

    os.chdir((MAIN_DIRECTORY + '/MERGnetcdfCEs/'))
    temporalRes = 3 #3 hours for TRMM

    #sort files
    files = filter(os.path.isfile, glob.glob('*.nc'))
    files.sort(key=lambda x: os.path.getmtime(x))

    for afile in files:
        fullFname = os.path.splitext(afile)[0]
        noFrameExtension = (fullFname.replace('_', '')).split('F')[0]
        ceUniqueID = 'F' +(fullFname.replace('_', '')).split('F')[1]
        fileDateTimeChar = (noFrameExtension.replace(':', '')).split('s')[1]
        fileDateTime = fileDateTimeChar.replace('-', '')
        fileDate = fileDateTime[:-6]
        fileHr1 = fileDateTime[-6:-4]

        cloudElementData = Dataset(afile, 'r', format='NETCDF4')
        brightnesstemp1 = cloudElementData.variables['brightnesstemp'][:,:,:]
        #latsrawCloudElements = cloudElementData.variables['latitude'][:]
        #lonsrawCloudElements = cloudElementData.variables['longitude'][:]

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
        precipRate = TRMMData.variables['pcp'][:,:,:]
        latsrawTRMMData = TRMMData.variables['latitude'][:]
        lonsrawTRMMData = TRMMData.variables['longitude'][:]
        lonsrawTRMMData[lonsrawTRMMData > 180] = lonsrawTRMMData[lonsrawTRMMData>180] - 360.
        LONTRMM, LATTRMM = np.meshgrid(lonsrawTRMMData, latsrawTRMMData)

        # TODO: From this point onward 'LAT' and 'LON' are used but undefined.
        #       perhaps they have been replaced with LONTRMM and LATTRMM??
        nygrd = len(LAT[:, 0])
        nxgrd = len(LON[0, :])

        precipRateMasked = ma.masked_array(precipRate, mask=(precipRate < 0.0))
        #---------regrid the TRMM data to the MERG dataset ----------------------------------
        #regrid using the do_regrid stuff from the Apache OCW
        regriddedTRMM = ma.zeros((0, nygrd, nxgrd))
        regriddedTRMM = utils.do_regrid(precipRateMasked[0,:,:], LATTRMM, LONTRMM, LAT, LON, order=1, mdi= -999999999)
        #----------------------------------------------------------------------------------

        TRMMData.close()


        # ------ NETCDF File stuff ------------------------------------
        thisFileName = MAIN_DIRECTORY+'/TRMMnetcdfCEs/'+ fileDateTime + ceUniqueID  + '.nc'
        currNetCDFTRMMData = Dataset(thisFileName, 'w', format='NETCDF4')
        currNetCDFTRMMData.description = 'Cloud Element ' + ceUniqueID + ' rainfall data'
        currNetCDFTRMMData.calendar = 'standard'
        currNetCDFTRMMData.conventions = 'COARDS'
        # dimensions
        currNetCDFTRMMData.createDimension('time', None)
        currNetCDFTRMMData.createDimension('lat', len(LAT[:,0]))
        currNetCDFTRMMData.createDimension('lon', len(LON[0,:]))
        # variables
        TRMMprecip = ('time', 'lat', 'lon',)
        times = currNetCDFTRMMData.createVariable('time', 'f8', ('time',))
        times.units = 'hours since '+ fileDateTime[:-6]
        latitude = currNetCDFTRMMData.createVariable('latitude', 'f8', ('lat',))
        longitude = currNetCDFTRMMData.createVariable('longitude', 'f8', ('lon',))
        rainFallacc = currNetCDFTRMMData.createVariable('precipitation_Accumulation', 'f8', TRMMprecip)
        rainFallacc.units = 'mm'

        longitude[:] = LON[0,:]
        longitude.units = 'degrees_east'
        longitude.long_name = 'Longitude'

        latitude[:] = LAT[:,0]
        latitude.units = 'degrees_north'
        latitude.long_name = 'Latitude'

        finalCETRMMvalues = ma.zeros((brightnesstemp1.shape))
        #-----------End most of NETCDF file stuff ------------------------------------
        for index, value in np.ndenumerate(brightnesstemp):
            latIndex, lonIndex = index
            if value > 0:
               finalCETRMMvalues[0,latIndex,lonIndex] = regriddedTRMM[int(np.where(LAT[:,0]==LAT[latIndex,0])[0]),\
                   int(np.where(LON[0,:]==LON[0,lonIndex])[0])]

        rainFallacc[:] = finalCETRMMvalues
        currNetCDFTRMMData.close()

        for index, value in np.ndenumerate(finalCETRMMvalues):
            precipTotal += value

        TRMMnumOfBoxes = np.count_nonzero(finalCETRMMvalues)
        TRMMArea = TRMMnumOfBoxes * XRES * YRES

        try:
            minCEprecipRate = np.min(finalCETRMMvalues[np.nonzero(finalCETRMMvalues)])
        except:
            minCEprecipRate = 0.0

        try:
            maxCEprecipRate = np.max(finalCETRMMvalues[np.nonzero(finalCETRMMvalues)])
        except:
            maxCEprecipRate = 0.0

        #add info to CLOUDELEMENTSGRAPH
        #TODO try block
        for eachdict in CLOUD_ELEMENT_GRAPH.nodes(ceUniqueID):
            if eachdict[1]['uniqueID'] == ceUniqueID:
                if not 'cloudElementPrecipTotal' in eachdict[1].keys():
                    eachdict[1]['cloudElementPrecipTotal'] = precipTotal
                if not 'cloudElementLatLonTRMM' in eachdict[1].keys():
                    eachdict[1]['cloudElementLatLonTRMM'] = finalCETRMMvalues
                if not 'TRMMArea' in eachdict[1].keys():
                    eachdict[1]['TRMMArea'] = TRMMArea
                if not 'CETRMMmin' in eachdict[1].keys():
                    eachdict[1]['CETRMMmin'] = minCEprecipRate
                if not 'CETRMMmax' in eachdict[1].keys():
                    eachdict[1]['CETRMMmax'] = maxCEprecipRate

        #clean up
        precipTotal = 0.0
        latsrawTRMMData = []
        lonsrawTRMMData = []
        #latsrawCloudElements=[]
        #lonsrawCloudElements=[]
        finalCETRMMvalues = []
        #CEprecipRate =[]
        brightnesstemp = []

    return allCEnodesTRMMdata
#**********************************************************************************************************************
def find_cloud_clusters(CEGraph):
    '''
    Purpose:: Determines the cloud clusters properties from the subgraphs in
        the graph i.e. prunes the graph according to the minimum depth

    Inputs:: CEGraph: a Networkx directed graph of the CEs with weighted edges
        according the area overlap between nodes (CEs) of consectuive frames

    Returns:: PRUNED_GRAPH: a Networkx directed graph of with CCs/ MCSs

    '''

    seenNode = []
    
    cloudClustersFile = open((MAIN_DIRECTORY + '/textFiles/cloudClusters.txt'), 'wb')

    for eachNode in CEGraph:
        #check if the node has been seen before
        if eachNode not in dict(enumerate(zip(*seenNode))):
            #look for all trees associated with node as the root
            thisPathDistanceAndLength = nx.single_source_dijkstra(CEGraph, eachNode)
            #determine the actual shortestPath and minimum depth/length
            maxDepthAndMinPath = find_max_depth_and_min_path(thisPathDistanceAndLength)
            if maxDepthAndMinPath:
                maxPathLength = maxDepthAndMinPath[0]
                shortestPath = maxDepthAndMinPath[1]

                #add nodes and paths to PRUNED_GRAPH
                for i in xrange(len(shortestPath)):
                    if PRUNED_GRAPH.has_node(shortestPath[i]) is False:
                        PRUNED_GRAPH.add_node(shortestPath[i])

                    #add edge if necessary
                    if i < (len(shortestPath)-1) and PRUNED_GRAPH.has_edge(shortestPath[i], shortestPath[i+1]) is False:
                        prunedGraphEdgeweight = CEGraph.get_edge_data(shortestPath[i], shortestPath[i+1])['weight']
                        PRUNED_GRAPH.add_edge(shortestPath[i], shortestPath[i+1], weight=prunedGraphEdgeweight)

                #note information in a file for consideration later i.e. checking to see if it works
                cloudClustersFile.write('\nSubtree pathlength is %d and path is %s' %(maxPathLength, shortestPath))
                #update seenNode info
                seenNode.append(shortestPath)

    print 'pruned graph'
    print 'number of nodes are: ', PRUNED_GRAPH.number_of_nodes()
    print 'number of edges are: ', PRUNED_GRAPH.number_of_edges()
    print ('*'*80)

    # graphTitle = 'Cloud Clusters observed over somewhere during sometime'
    # plotting.draw_graph(PRUNED_GRAPH, graphTitle, MAIN_DIRECTORY, edgeWeight)
    cloudClustersFile.close()

    return PRUNED_GRAPH
#**********************************************************************************************************************
def find_MCC(prunedGraph):
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

    #maxShieldNode = ''
    orderedPath = []
    treeTraversalList = []
    definiteMCCFlag = False
    unDirGraph = nx.Graph()
    aSubGraph = nx.DiGraph()
    #definiteMCSFlag = False


    #connected_components is not available for DiGraph, so generate graph as undirected
    unDirGraph = PRUNED_GRAPH.to_undirected()
    subGraph = nx.connected_component_subgraphs(unDirGraph)

    #for each path in the subgraphs determined
    for path in subGraph:
        #definite is a subTree provided the duration is longer than 3 hours
        if len(path.nodes()) > MIN_MCS_DURATION:
            orderedPath = path.nodes()
            orderedPath.sort(key=lambda item: (len(item.split('C')[0]), item.split('C')[0]))
            #definiteMCS.append(orderedPath)

            #build back DiGraph for checking purposes/paper purposes
            aSubGraph.add_nodes_from(path.nodes())
            for eachNode in path.nodes():
                if prunedGraph.predecessors(eachNode):
                    for node in prunedGraph.predecessors(eachNode):
                        aSubGraph.add_edge(node, eachNode, weight=edgeWeight[0])

                if prunedGraph.successors(eachNode):
                    for node in prunedGraph.successors(eachNode):
                        aSubGraph.add_edge(eachNode, node, weight=edgeWeight[0])
            imgTitle = 'CC'+str(imgCount+1)
            # plotting.draw_graph(aSubGraph, imgTitle, MAIN_DIRECTORY, edgeWeight) #for eachNode in path:
            imgCount += 1
            #----------end build back ---------------------------------------------
            mergeList, splitList = has_merges_or_splits(path)
            #add node behavior regarding neutral, merge, split or both
            for node in path:
                if node in mergeList and node in splitList:
                    add_node_behavior_identifier(node, 'B')
                elif node in mergeList and not node in splitList:
                    add_node_behavior_identifier(node, 'M')
                elif node in splitList and not node in mergeList:
                    add_node_behavior_identifier(node, 'S')
                else:
                    add_node_behavior_identifier(node, 'N')
            #Do the first part of checking for the MCC feature
            #find the path
            treeTraversalList = traverse_tree(aSubGraph, orderedPath[0], [], [])
            #print 'treeTraversalList is ', treeTraversalList
            #check the nodes to determine if a MCC on just the area criteria (consecutive nodes meeting the area and temp requirements)
            MCCList = checked_nodes_MCC(prunedGraph, treeTraversalList)
            for aDict in MCCList:
                for eachNode in aDict['fullMCSMCC']:
                    add_node_MCS_identifier(eachNode[0], eachNode[1])

            #do check for if MCCs overlap
            if MCCList:
                if len(MCCList) > 1:
                    for count in range(len(MCCList)): #for eachDict in MCCList:
                        #if there are more than two lists
                        if count >= 1:
                            #and the first node in this list
                            eachList = list(x[0] for x in MCCList[count]['possMCCList'])
                            eachList.sort(key=lambda nodeID: (len(nodeID.split('C')[0]), nodeID.split('C')[0]))
                            if eachList:
                                fNode = eachList[0]
                                #get the lastNode in the previous possMCC list
                                eachList = list(x[0] for x in MCCList[(count-1)]['possMCCList'])
                                eachList.sort(key=lambda nodeID: (len(nodeID.split('C')[0]), nodeID.split('C')[0]))
                                if eachList:
                                    lNode = eachList[-1]
                                    if lNode in CLOUD_ELEMENT_GRAPH.predecessors(fNode):
                                        for aNode in CLOUD_ELEMENT_GRAPH.predecessors(fNode):
                                            if aNode in eachList and aNode == lNode:
                                            #if edge_data is equal or less than to the exisitng edge in the tree append one to the other
                                                if CLOUD_ELEMENT_GRAPH.get_edge_data(aNode,fNode)['weight'] <= \
                                                        CLOUD_ELEMENT_GRAPH.get_edge_data(lNode,fNode)['weight']:
                                                    MCCList[count-1]['possMCCList'].extend(MCCList[count]['possMCCList'])
                                                    MCCList[count-1]['fullMCSMCC'].extend(MCCList[count]['fullMCSMCC'])
                                                    MCCList[count-1]['durationAandB'] += MCCList[count]['durationAandB']
                                                    MCCList[count-1]['CounterCriteriaA'] += MCCList[count]['CounterCriteriaA']
                                                    MCCList[count-1]['highestMCCnode'] = MCCList[count]['highestMCCnode']
                                                    MCCList[count-1]['frameNum'] = MCCList[count]['frameNum']
                                                    removeList.append(count)
                #update the MCCList
                if removeList:
                    for i in removeList:
                        if (len(MCCList) - 1) > i:
                            del MCCList[i]
                            removeList = []

            #check if the nodes also meet the duration criteria and the shape crieria
            for eachDict in MCCList:
                #order the fullMCSMCC list, then run maximum extent and eccentricity criteria
                if (eachDict['durationAandB'] * TRES) >= MINIMUM_DURATION and \
                        (eachDict['durationAandB'] * TRES) <= MAXIMUM_DURATION:
                    eachList = list(x[0] for x in eachDict['fullMCSMCC'])
                    eachList.sort(key=lambda nodeID: (len(nodeID.split('C')[0]), nodeID.split('C')[0]))
                    eachMCCList = list(x[0] for x in eachDict['possMCCList'])
                    eachMCCList.sort(key=lambda nodeID: (len(nodeID.split('C')[0]), nodeID.split('C')[0]))

                    #update the nodemcsidentifer behavior
                    #find the first element eachMCCList in eachList, and ensure everything ahead of it is indicated as 'I',
                    #find last element in eachMCCList in eachList and ensure everything after it is indicated as 'D'
                    #ensure that everything between is listed as 'M'
                    for eachNode in eachList[:(eachList.index(eachMCCList[0]))]:
                        add_node_MCS_identifier(eachNode, 'I')

                    add_node_MCS_identifier(eachMCCList[0], 'M')

                    for eachNode in eachList[(eachList.index(eachMCCList[-1])+1):]:
                        add_node_MCS_identifier(eachNode, 'D')

                    #update definiteMCS list
                    for eachNode in orderedPath[(orderedPath.index(eachMCCList[-1])+1):]:
                        add_node_MCS_identifier(eachNode, 'D')

                    #run maximum extent and eccentricity criteria
                    _, definiteMCCFlag = max_extent_and_eccentricity(eachList)
                    #maxExtentNode, definiteMCCFlag = max_extent_and_eccentricity(eachList)
                    #print 'maxExtentNode, definiteMCCFlag ', maxExtentNode, definiteMCCFlag
                    if definiteMCCFlag == True:
                        definiteMCC.append(eachList)


            definiteMCS.append(orderedPath)

            #reset for next subGraph
            aSubGraph.clear()
            orderedPath = []
            MCCList = []
            #MCSList =[]
            #definiteMCSFlag = False

    return definiteMCC, definiteMCS
#**********************************************************************************************************************
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

    #check one level infront first...if something does exisit, stick it at the front of the stack
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
#**********************************************************************************************************************
def checked_nodes_MCC(prunedGraph, nodeList):
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
    thisdict = {} #will have the same items as the cloudElementDict
    cloudElementAreaB = 0.0
    #cloudElementAreaA = 0.0
    #epsilon = 0.0
    #frameNum =0
    oldNode = ''
    potentialMCCList = []
    #durationAandB = 0

    #check for if the list contains only one string/node
    if type(nodeList) is str:
        oldNode = nodeList
        nodeList = []
        nodeList.append(oldNode)

    for node in nodeList:
        thisdict = this_dict(node)
        counterCriteriaAFlag = False
        counterCriteriaBFlag = False
        #existingFrameFlag = False

        if thisdict['cloudElementArea'] >= OUTER_CLOUD_SHIELD_AREA:
            counterCriteriaAFlag = True
            initiationFlag = True
            maturityFlag = False

            #check if criteriaA is met
            cloudElementAreaA, _ = check_criteria(thisdict['cloudElementLatLon'], OUTER_CLOUD_SHIELD_TEMPERATURE, False)
            #cloudElementAreaA, criteriaA = check_criteria(thisdict['cloudElementLatLon'], OUTER_CLOUD_SHIELD_TEMPERATURE)
            #TODO: calcuate the eccentricity at this point and read over????or create a new field in the dict

            if cloudElementAreaA >= OUTER_CLOUD_SHIELD_AREA:
                #check if criteriaB is met
                cloudElementAreaB, criteriaB = check_criteria(thisdict['cloudElementLatLon'], INNER_CLOUD_SHIELD_TEMPERATURE, True)

                #if Criteria A and B have been met, then the MCC is initiated, i.e. store node as potentialMCC
                if cloudElementAreaB >= INNER_CLOUD_SHIELD_AREA:
                    #TODO: add another field to the dictionary for the OUTER_AREA_SHIELD area
                    counterCriteriaBFlag = True
                    #append this information on to the dictionary
                    add_info_this_dict(node, cloudElementAreaB, criteriaB)
                    initiationFlag = False
                    maturityFlag = True
                    stage = 'M'
                    potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,\
                                         counterCriteriaAFlag, counterCriteriaBFlag)
                else:
                    #criteria B failed
                    counterCriteriaBFlag = False
                    if initiationFlag == True:
                        stage = 'I'
                        potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,\
                                             counterCriteriaAFlag, counterCriteriaBFlag)

                    elif (initiationFlag == False and maturityFlag == True) or decayFlag == True:
                        decayFlag = True
                        maturityFlag = False
                        stage = 'D'
                        potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,\
                                             counterCriteriaAFlag, counterCriteriaBFlag)
            else:
                #criteria A failed
                counterCriteriaAFlag = False
                counterCriteriaBFlag = False
                #add as a CE before or after the main feature
                if initiationFlag == True or (initiationFlag == False and maturityFlag == True):
                    stage = 'I'
                    potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,\
                                         counterCriteriaAFlag, counterCriteriaBFlag)
                elif (initiationFlag == False and maturityFlag == False) or decayFlag == True:
                    stage = 'D'
                    decayFlag = True
                    potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,\
                                         counterCriteriaAFlag, counterCriteriaBFlag)
                elif (initiationFlag == False and maturityFlag == False and decayFlag == False):
                    stage = 'I'
                    potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,\
                                         counterCriteriaAFlag, counterCriteriaBFlag)
        else:
            #criteria A failed
            counterCriteriaAFlag = False
            counterCriteriaBFlag = False
            #add as a CE before or after the main feature
            if initiationFlag == True or (initiationFlag == False and maturityFlag == True):
                stage = 'I'
                potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,\
                                     counterCriteriaAFlag, counterCriteriaBFlag)
            elif (initiationFlag == False and maturityFlag == False) or decayFlag == True:
                stage = 'D'
                decayFlag = True
                potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,\
                                     counterCriteriaAFlag, counterCriteriaBFlag)
            elif (initiationFlag == False and maturityFlag == False and decayFlag == False):
                stage = 'I'
                potentialMCCList = update_MCC_list(prunedGraph, potentialMCCList, node, stage,\
                                     counterCriteriaAFlag, counterCriteriaBFlag)
    return potentialMCCList
#**********************************************************************************************************************
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
        #list empty
        stage = 'I'
        if counterCriteriaAFlag == True and counterCriteriaBFlag == True:
            potentialMCCList.append({'possMCCList':[(node,stage)], 'fullMCSMCC':[(node,stage)], 'CounterCriteriaA': 1,\
                 'durationAandB': 1, 'highestMCCnode':node, 'frameNum':frameNum})
        elif counterCriteriaAFlag == True and counterCriteriaBFlag == False:
            potentialMCCList.append({'possMCCList':[], 'fullMCSMCC':[(node, stage)], 'CounterCriteriaA': 1,\
                 'durationAandB': 0, 'highestMCCnode':'', 'frameNum':0})
        elif counterCriteriaAFlag == False and counterCriteriaBFlag == False:
            potentialMCCList.append({'possMCCList':[], 'fullMCSMCC':[(node, stage)], 'CounterCriteriaA': 0,\
                 'durationAandB': 0, 'highestMCCnode':'', 'frameNum':0})
    else:
        #list not empty
        predecessorsFlag, index = is_there_a_link(prunedGraph, 1, node, potentialMCCList, 1)
        if predecessorsFlag == True:
            for eachNode in potentialMCCList[index]['possMCCList']:
                if int((eachNode[0].split('CE')[0]).split('F')[1]) == frameNum :
                    existingFrameFlag = True
            #this MUST come after the check for the existing frame
            if counterCriteriaAFlag == True and counterCriteriaBFlag == True:
                stage = 'M'
                potentialMCCList[index]['possMCCList'].append((node, stage))
                potentialMCCList[index]['fullMCSMCC'].append((node, stage))


            if existingFrameFlag == False:
                if counterCriteriaAFlag == True and counterCriteriaBFlag == True:
                    stage = 'M'
                    potentialMCCList[index]['CounterCriteriaA'] += 1
                    potentialMCCList[index]['durationAandB'] += 1
                    if frameNum > potentialMCCList[index]['frameNum']:
                        potentialMCCList[index]['frameNum'] = frameNum
                        potentialMCCList[index]['highestMCCnode'] = node
                    return potentialMCCList

                #if this frameNum doesn't exist and this frameNum is less than the MCC node max frame Num (including 0), then append to fullMCSMCC list
                if frameNum > potentialMCCList[index]['frameNum'] or potentialMCCList[index]['frameNum'] == 0:
                    stage = 'I'
                    if counterCriteriaAFlag == True and counterCriteriaBFlag == False:
                        potentialMCCList.append({'possMCCList':[], 'fullMCSMCC':[(node, stage)], 'CounterCriteriaA': 1,\
                             'durationAandB': 0, 'highestMCCnode':'', 'frameNum':0})
                        return potentialMCCList
                    elif counterCriteriaAFlag == False and counterCriteriaBFlag == False:
                        potentialMCCList.append({'possMCCList':[], 'fullMCSMCC':[(node, stage)], 'CounterCriteriaA': 0,\
                             'durationAandB': 0, 'highestMCCnode':'', 'frameNum':0})
                        return potentialMCCList

            #if predecessor and this frame number already exist in the MCC list, add the current node to the fullMCSMCC list
            if existingFrameFlag == True:
                if counterCriteriaAFlag == True and counterCriteriaBFlag == False:
                    potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                    potentialMCCList[index]['CounterCriteriaA'] += 1
                    return potentialMCCList
                if counterCriteriaAFlag == False:
                    potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                    return potentialMCCList

        if predecessorsFlag == False:
            successorsFlag, index = is_there_a_link(prunedGraph, 2, node, potentialMCCList, 2)

            if successorsFlag == True:
                for eachNode in potentialMCCList[index]['possMCCList']:
                    if int((eachNode[0].split('CE')[0]).split('F')[1]) == frameNum:
                        existingFrameFlag = True

                if counterCriteriaAFlag == True and counterCriteriaBFlag == True:
                    stage = 'M'
                    potentialMCCList[index]['possMCCList'].append((node, stage))
                    potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                    if frameNum > potentialMCCList[index]['frameNum'] or potentialMCCList[index]['frameNum'] == 0:
                        potentialMCCList[index]['frameNum'] = frameNum
                        potentialMCCList[index]['highestMCCnode'] = node
                    return potentialMCCList


                if existingFrameFlag == False:
                    if stage == 'M':
                        stage = 'D'
                    if counterCriteriaAFlag == True and counterCriteriaBFlag == True:
                        potentialMCCList[index]['CounterCriteriaA'] += 1
                        potentialMCCList[index]['durationAandB'] += 1
                    elif counterCriteriaAFlag == True:
                        potentialMCCList[index]['CounterCriteriaA'] += 1
                    elif counterCriteriaAFlag == False:
                        potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                        return potentialMCCList
    #if predecessor and this frame number already exist in the MCC list, add the current node to the fullMCSMCC list
                else:
                    if counterCriteriaAFlag == True and counterCriteriaBFlag == False:
                        potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                        potentialMCCList[index]['CounterCriteriaA'] += 1
                        return potentialMCCList
                    if counterCriteriaAFlag == False:
                        potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                        return potentialMCCList

        #if this node isn't connected to exisiting MCCs check if it is connected to exisiting MCSs ...
        if predecessorsFlag == False and successorsFlag == False:
            stage = 'I'
            predecessorsMCSFlag, index = is_there_a_link(prunedGraph, 1, node, potentialMCCList, 2)
            if predecessorsMCSFlag == True:
                if counterCriteriaAFlag == True and counterCriteriaBFlag == True:
                    potentialMCCList[index]['possMCCList'].append((node, 'M'))
                    potentialMCCList[index]['fullMCSMCC'].append((node, 'M'))
                    potentialMCCList[index]['durationAandB'] += 1
                    if frameNum > potentialMCCList[index]['frameNum']:
                        potentialMCCList[index]['frameNum'] = frameNum
                        potentialMCCList[index]['highestMCCnode'] = node
                    return potentialMCCList

                if potentialMCCList[index]['frameNum'] == 0 or frameNum <= potentialMCCList[index]['frameNum']:
                    if counterCriteriaAFlag == True and counterCriteriaBFlag == False:
                        potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                        potentialMCCList[index]['CounterCriteriaA'] += 1
                        return potentialMCCList
                    elif counterCriteriaAFlag == False:
                        potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                        return potentialMCCList
            else:
                successorsMCSFlag, index = is_there_a_link(prunedGraph, 2, node, potentialMCCList, 2)
                if successorsMCSFlag == True:
                    if counterCriteriaAFlag == True and counterCriteriaBFlag == True:
                        potentialMCCList[index]['possMCCList'].append((node, 'M'))
                        potentialMCCList[index]['fullMCSMCC'].append((node, 'M'))
                        potentialMCCList[index]['durationAandB'] += 1
                        if frameNum > potentialMCCList[index]['frameNum']:
                            potentialMCCList[index]['frameNum'] = frameNum
                            potentialMCCList[index]['highestMCCnode'] = node
                        return potentialMCCList


                    if potentialMCCList[index]['frameNum'] == 0 or frameNum <= potentialMCCList[index]['frameNum']:
                        if counterCriteriaAFlag == True and counterCriteriaBFlag == False:
                            potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                            potentialMCCList[index]['CounterCriteriaA'] += 1
                            return potentialMCCList
                        elif counterCriteriaAFlag == False:
                            potentialMCCList[index]['fullMCSMCC'].append((node, stage))
                            return potentialMCCList

            #if this node isn't connected to existing MCCs or MCSs, create a new one ...
            if predecessorsFlag == False and predecessorsMCSFlag == False and successorsFlag == False and\
                     successorsMCSFlag == False:
                if counterCriteriaAFlag == True and counterCriteriaBFlag == True:
                    potentialMCCList.append({'possMCCList':[(node, stage)], 'fullMCSMCC':[(node, stage)],\
                         'CounterCriteriaA': 1, 'durationAandB': 1, 'highestMCCnode':node, 'frameNum':frameNum})
                elif counterCriteriaAFlag == True and counterCriteriaBFlag == False:
                    potentialMCCList.append({'possMCCList':[], 'fullMCSMCC':[(node, stage)], 'CounterCriteriaA': 1,\
                         'durationAandB': 0, 'highestMCCnode':'', 'frameNum':0})
                elif counterCriteriaAFlag == False and counterCriteriaBFlag == False:
                    potentialMCCList.append({'possMCCList':[], 'fullMCSMCC':[(node, stage)], 'CounterCriteriaA': 0,\
                         'durationAandB': 0, 'highestMCCnode':'', 'frameNum':0})

    return potentialMCCList
#**********************************************************************************************************************
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

    #check parents
    if upOrDown == 1:
        for aNode in prunedGraph.predecessors(node):
            #reset the index counter for this node search through potentialMCCList
            index = -1
            for mccDict in potentialMCCList:
                index += 1
                if aNode in list(x[0] for x in mccDict[checkList]):
                    thisFlag = True
                    #get out of looping so as to avoid the flag being written over when another node in the predecesor 
                    #list is checked
                    return thisFlag, index

    #check children
    if upOrDown == 2:
        for aNode in prunedGraph.successors(node):
            #reset the index counter for this node search through potentialMCCList
            index = -1
            for mccDict in potentialMCCList:
                index += 1

                if aNode in list(x[0] for x in mccDict[checkList]):
                    thisFlag = True
                    return thisFlag, index

    return thisFlag, index
#**********************************************************************************************************************
def max_extent_and_eccentricity(eachList):
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
    #maxShieldEccentricity = 0.0
    definiteMCCFlag = False

    if eachList:
        for eachNode in eachList:
            if (this_dict(eachNode)['nodeMCSIdentifier'] == 'M' or this_dict(eachNode)['nodeMCSIdentifier'] == 'D') and\
                    this_dict(eachNode)['cloudElementArea'] > maxShieldArea:
                maxShieldNode = eachNode
                maxShieldArea = this_dict(eachNode)['cloudElementArea']

        #maxShieldEccentricity = this_dict(maxShieldNode)['cloudElementEccentricity']
        if this_dict(maxShieldNode)['cloudElementEccentricity'] >= ECCENTRICITY_THRESHOLD_MIN and \
                this_dict(maxShieldNode)['cloudElementEccentricity'] <= ECCENTRICITY_THRESHOLD_MAX:
            #criteria met
            definiteMCCFlag = True

    return maxShieldNode, definiteMCCFlag
#**********************************************************************************************************************
def find_max_depth_and_min_path(thisPathDistanceAndLength):
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

    #maxPathLength for the node in question
    maxPathLength = max(len (values) for values in thisPathDistanceAndLength[1].values())

    #if the duration is shorter then the min MCS length, then don't store!
    if maxPathLength < MIN_MCS_DURATION: #MINIMUM_DURATION :
        minDistanceAndMaxPath = ()

    #else find the min path and max depth
    else:
        #max path distance for the node in question
        minPath = max(values for values in thisPathDistanceAndLength[0].values())

        #check to determine the shortest path from the longest paths returned
        for pathDistance, path in itertools.izip(thisPathDistanceAndLength[0].values(), thisPathDistanceAndLength[1].values()):
            pathLength = len(path)
            #if pathLength is the same as the maxPathLength, then look the pathDistance to determine if the min
            if pathLength == maxPathLength:
                if pathDistance <= minPath:
                    minPath = pathLength
                    #store details if absolute minPath and deepest
                    minDistanceAndMaxPath = (pathDistance, path)
    return minDistanceAndMaxPath
#**********************************************************************************************************************
def this_dict(thisNode):
    '''
    Purpose::
        Return dictionary from graph if node exist in tree

    Input::
        thisNode: a string representing the CE to get the information for

    Output ::
        eachdict[1]: a dictionary representing the info associated with thisNode from the graph

    '''
    for eachdict in CLOUD_ELEMENT_GRAPH.nodes(thisNode):
        if eachdict[1]['uniqueID'] == thisNode:
            return eachdict[1]
#**********************************************************************************************************************
def check_criteria(thisCloudElementLatLon, aTemperature, bFlag):
    '''
    Purpose:: Determine if criteria B is met for a CEGraph

    Inputs:: thisCloudElementLatLon: 2D array of (lat,lon) variable from the node dictionary being currently considered
        aTemperature:a integer representing the temperature maximum for masking
        bFlag: a boolean indicating if we need to care about cloudElementCriteriaBLatLon

    Returns:: cloudElementArea: a floating-point number representing the area in the array that meet the criteria 

    '''
    cloudElementCriteriaBLatLon = []

    _, ceCounter = ndimage.measurements.label(thisCloudElementLatLon, structure=STRUCTURING_ELEMENT)
    #frame, ceCounter = ndimage.measurements.label(thisCloudElementLatLon, structure=STRUCTURING_ELEMENT)
    #frameceCounter = 0
    #determine min and max values in lat and lon, then use this to generate teh array from LAT,LON meshgrid

    minLat = min(x[0] for x in thisCloudElementLatLon)
    maxLat = max(x[0] for x in thisCloudElementLatLon)
    minLon = min(x[1] for x in thisCloudElementLatLon)
    maxLon = max(x[1] for x in thisCloudElementLatLon)

    minLatIndex = np.argmax(LAT[:,0] == minLat)
    maxLatIndex = np.argmax(LAT[:,0] == maxLat)
    minLonIndex = np.argmax(LON[0,:] == minLon)
    maxLonIndex = np.argmax(LON[0,:] == maxLon)

    criteriaBframe = ma.zeros(((abs(maxLatIndex - minLatIndex)+1), (abs(maxLonIndex - minLonIndex)+1)))

    for x in thisCloudElementLatLon:
        #to store the values of the subset in the new array, remove the minLatIndex and minLonindex from the
        #index given in the original array to get the indices for the new array
        criteriaBframe[(np.argmax(LAT[:,0] == x[0]) - minLatIndex), (np.argmax(LON[0,:] == x[1]) - minLonIndex)] = x[2]

    #keep only those values < aTemperature
    tempMask = ma.masked_array(criteriaBframe, mask=(criteriaBframe >= aTemperature), fill_value = 0)

    #get the actual values that the mask returned
    criteriaB = ma.zeros((criteriaBframe.shape)).astype('int16')

    for index, value in utils.maenumerate(tempMask):
        latIndex, lonIndex = index
        criteriaB[latIndex, lonIndex] = value

    for _ in xrange(ceCounter):
        #[0] is time dimension. Determine the actual values from the data
        #loc is a masked array
        #***** returns elements down then across thus (6,4) is 6 arrays deep of size 4
        try:

            loc = ndimage.find_objects(criteriaB)[0]
        except:
            #this would mean that no objects were found meeting criteria B
            print 'no objects at this temperature!'
            cloudElementArea = 0.0
            return cloudElementArea, cloudElementCriteriaBLatLon
        try:
            cloudElementCriteriaB = ma.zeros((criteriaB.shape))
            cloudElementCriteriaB = criteriaB[loc]
        except:
            print 'YIKESS'
            print 'ceCounter ', ceCounter, criteriaB.shape
            print 'criteriaB ', criteriaB

        if bFlag == True:
            for index, value in np.ndenumerate(cloudElementCriteriaB):
                if value != 0:
                    _, lat, lon = index
                    #t,lat,lon = index
                    #add back on the minLatIndex and minLonIndex to find the true lat, lon values
                    latLonTuple = (LAT[(lat),0], LON[0,(lon)], value)
                    cloudElementCriteriaBLatLon.append(latLonTuple)

        cloudElementArea = np.count_nonzero(cloudElementCriteriaB) * XRES * YRES
        #do some cleaning up
        tempMask = []
        criteriaB = []
        cloudElementCriteriaB = []

        return cloudElementArea, cloudElementCriteriaBLatLon
#**********************************************************************************************************************
def has_merges_or_splits(nodeList):
    '''
    Purpose:: Determine if nodes within a path defined from shortest_path splittingNodeDict
    Inputs:: nodeList: list of strings representing the nodes from a path
    Returns:: splitList: a list of strings representing all the nodes in the path that split
        mergeList: a list of strings representing all the nodes in the path that merged
    '''
    mergeList = []
    splitList = []

    for node,numParents in PRUNED_GRAPH.in_degree(nodeList).items():
        if numParents > 1:
            mergeList.append(node)

    for node, numChildren in PRUNED_GRAPH.out_degree(nodeList).items():
        if numChildren > 1:
            splitList.append(node)
    #sort
    splitList.sort(key=lambda item: (len(item.split('C')[0]), item.split('C')[0]))
    mergeList.sort(key=lambda item: (len(item.split('C')[0]), item.split('C')[0]))

    return mergeList, splitList
#**********************************************************************************************************************
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
#**********************************************************************************************************************
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
            #i.e. PRUNED_GRAPH.predecessors(aNode) is empty
            return path, numOfChildren
    except:
        #i.e. PRUNED_GRAPH.predecessors(aNode) threw an exception
        return path, numOfChildren
#**********************************************************************************************************************
def add_info_this_dict(thisNode, cloudElementArea, criteriaB):
    '''
    Purpose:: Update original dictionary node with information

    Inputs:: thisNode: a string representing the unique ID of a node
        cloudElementArea: a floating-point number representing the area of the cloud element
        criteriaB: a masked array of floating-point numbers representing the lat,lons meeting the criteria

    Returns:: None

    '''
    for eachdict in CLOUD_ELEMENT_GRAPH.nodes(thisNode):
        if eachdict[1]['uniqueID'] == thisNode:
            eachdict[1]['CriteriaBArea'] = cloudElementArea
            eachdict[1]['CriteriaBLatLon'] = criteriaB
    return
#**********************************************************************************************************************
def add_node_behavior_identifier(thisNode, nodeBehaviorIdentifier):
    '''
    Purpose:: add an identifier to the node dictionary to indicate splitting, merging or neither node

    Inputs:: thisNode: a string representing the unique ID of a node
        nodeBehaviorIdentifier: a string representing the behavior S- split, M- merge, B- both split and merge, 
        N- neither split or merge

    Returns:: None

    '''
    for eachdict in CLOUD_ELEMENT_GRAPH.nodes(thisNode):
        if eachdict[1]['uniqueID'] == thisNode:
            if not 'nodeBehaviorIdentifier' in eachdict[1].keys():
                eachdict[1]['nodeBehaviorIdentifier'] = nodeBehaviorIdentifier
    return
#**********************************************************************************************************************
def add_node_MCS_identifier(thisNode, nodeMCSIdentifier):
    '''
    Purpose:: Add an identifier to the node dictionary to indicate splitting, merging or neither node

    Inputs:: thisNode: a string representing the unique ID of a node
        nodeMCSIdentifier: a string representing the stage of the MCS lifecyle  'I' for Initiation, 'M' for Maturity,
        'D' for Decay

    Returns:: None

    '''
    for eachdict in CLOUD_ELEMENT_GRAPH.nodes(thisNode):
        if eachdict[1]['uniqueID'] == thisNode:
            if not 'nodeMCSIdentifier' in eachdict[1].keys():
                eachdict[1]['nodeMCSIdentifier'] = nodeMCSIdentifier
    return
#**********************************************************************************************************************
def update_node_MCS_identifier(thisNode, nodeMCSIdentifier):
    '''
    Purpose:: Update an identifier to the node dictionary to indicate splitting, merging or neither node

    Inputs:: thisNode: thisNode: a string representing the unique ID of a node
        nodeMCSIdentifier: a string representing the stage of the MCS lifecyle  'I' for Initiation, 'M' for Maturity, 'D' for Decay

    Returns:: None

    '''
    for eachdict in CLOUD_ELEMENT_GRAPH.nodes(thisNode):
        if eachdict[1]['uniqueID'] == thisNode:
            eachdict[1]['nodeMCSIdentifier'] = nodeMCSIdentifier

    return
#**********************************************************************************************************************
def eccentricity(cloudElementLatLon):
    '''
    Purpose:: Determines the eccentricity (shape) of contiguous boxes
        Values tending to 1 are more circular by definition, whereas
        values tending to 0 are more linear

    Inputs:: cloudElementLatLon: 2D array in (lat,lon) representing T_bb contiguous squares

    Returns:: epsilon: a floating-point representing the eccentricity of the matrix passed

    '''

    epsilon = 0.0

    #loop over all lons and determine longest (non-zero) col
    #loop over all lats and determine longest (non-zero) row

    for _ in cloudElementLatLon:
        #assign a matrix to determine the legit values

        nonEmptyLons = sum(sum(cloudElementLatLon) > 0)
        nonEmptyLats = sum(sum(cloudElementLatLon.transpose()) > 0)

        lonEigenvalues = 1.0 * nonEmptyLats / (nonEmptyLons + 0.001) #for long oval on y axis
        latEigenvalues = 1.0 * nonEmptyLons / (nonEmptyLats + 0.001) #for long oval on x-axs
        epsilon = min(latEigenvalues, lonEigenvalues)

    return epsilon
#**********************************************************************************************************************
def cloud_element_overlap(currentCELatLons, previousCELatLons):
    '''
    Purpose::
        Determines the percentage overlap between two list of lat-lons passed

    Input::
        currentCELatLons: a list of tuples for the current CE
        previousCELatLons: a list of tuples for the other CE being considered

    Output::
        percentageOverlap: a floating-point representing the number of overlapping lat_lon tuples
        areaOverlap: a floating-point number representing the area overlapping

    '''

    latlonprev = []
    latloncurr = []
    count = 0
    percentageOverlap = 0.0
    areaOverlap = 0.0

    #remove the temperature from the tuples for currentCELatLons and previousCELatLons then check for overlap
    latlonprev = [(x[0], x[1]) for x in previousCELatLons]
    latloncurr = [(x[0], x[1]) for x in currentCELatLons]

    #find overlap
    count = len(list(set(latloncurr) & set(latlonprev)))

    #find area overlap
    areaOverlap = count * XRES * YRES

    #find percentage
    percentageOverlap = max(((count * 1.0) / (len(latloncurr) * 1.0)), ((count * 1.0) / (len(latlonprev) * 1.0)))

    return percentageOverlap, areaOverlap
#**********************************************************************************************************************
