'''
# All the functions for the MCC search algorithm
# Following RCMES dataformat in format (t,lat,lon), value
'''

from datetime import timedelta
import glob
import itertools
from netCDF4 import Dataset, date2num
import numpy as np
import numpy.ma as ma
import os
from scipy import ndimage

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
XRES = 4.0				#x direction spatial resolution in km
YRES = 4.0				#y direction spatial resolution in km
TRES = 1 				#temporal resolution in hrs
LAT_DISTANCE = 111.0 	#the avg distance in km for 1deg lat for the region being considered 
LON_DISTANCE = 111.0    #the avg distance in km for 1deg lon for the region being considered
STRUCTURING_ELEMENT = [[0,1,0],[1,1,1],[0,1,0]] #the matrix for determining the pattern for the contiguous boxes and must
    											#have same rank of the matrix it is being compared against 
#criteria for determining cloud elements and edges
T_BB_MAX = 243  #warmest temp to allow (-30C to -55C according to Morel and Sensi 2002)
T_BB_MIN = 218  #cooler temp for the center of the system
CONVECTIVE_FRACTION = 0.90 #the min temp/max temp that would be expected in a CE.. this is highly conservative (only a 10K difference)
MIN_MCS_DURATION = 3    #minimum time for a MCS to exist
AREA_MIN = 2400.0		#minimum area for CE criteria in km^2 according to Vila et al. (2008) is 2400
MIN_OVERLAP= 10000.00   #km^2  from Williams and Houze 1987, indir ref in Arnaud et al 1992

#---the MCC criteria
ECCENTRICITY_THRESHOLD_MAX = 1.0  #tending to 1 is a circle e.g. hurricane, 
ECCENTRICITY_THRESHOLD_MIN = 0.70 #tending to 0 is a linear e.g. squall line
OUTER_CLOUD_SHIELD_AREA = 80000.0 #km^2
INNER_CLOUD_SHIELD_AREA = 30000.0 #km^2
OUTER_CLOUD_SHIELD_TEMPERATURE = 233 #in K
INNER_CLOUD_SHIELD_TEMPERATURE = 213 #in K
MINIMUM_DURATION = 6 #min number of frames the MCC must exist for (assuming hrly frames, MCCs is 6hrs)
MAXIMUM_DURATION = 24#max number of framce the MCC can last for 
#------------------- End user defined Variables -------------------
edgeWeight = [1,2,3] #weights for the graph edges
#graph object fo the CEs meeting the criteria
CLOUD_ELEMENT_GRAPH = nx.DiGraph()
#graph meeting the CC criteria
PRUNED_GRAPH = nx.DiGraph()
#get lat lons from iomethods.py

#------------------------ End GLOBAL VARS -------------------------
#************************ Begin Functions *************************
#******************************************************************
def findCloudElements(mergImgs,timelist,mainStrDir,LAT,LON,TRMMdirName=None):
	'''
	Purpose::
	    Determines the contiguous boxes for a given time of the satellite images i.e. each frame
        using scipy ndimage package
	
	Input::	
		mergImgs: masked numpy array in (time,lat,lon),T_bb representing the satellite data. This is masked based on the
		maximum acceptable temperature, T_BB_MAX
		timelist: a list of python datatimes
		LAT - a 2D array of the set of latitudes from the files opened
		LON - a 2D array of the set of longitudes from the files opened
		TRMMdirName (optional): string representing the path where to find the TRMM datafiles
		
	Output::
	    CLOUD_ELEMENT_GRAPH: a Networkx directed graph where each node contains the information in cloudElementDict
	    The nodes are determined according to the area of contiguous squares. The nodes are linked through weighted edges.

		cloudElementDict = {'uniqueID': unique tag for this CE, 
							'cloudElementTime': time of the CE,
							'cloudElementLatLon': (lat,lon,value) of MERG data of CE, 
							'cloudElementCenter':list of floating-point [lat,lon] representing the CE's center 
							'cloudElementArea':floating-point representing the area of the CE, 
							'cloudElementEccentricity': floating-point representing the shape of the CE, 
							'cloudElementTmax':integer representing the maximum Tb in CE, 
							'cloudElementTmin': integer representing the minimum Tb in CE, 
							'cloudElementPrecipTotal':floating-point representing the sum of all rainfall in CE if TRMMdirName entered,
							'cloudElementLatLonTRMM':(lat,lon,value) of TRMM data in CE if TRMMdirName entered, 
							'TRMMArea': floating-point representing the CE if TRMMdirName entered,
							'CETRMMmax':floating-point representing the max rate in the CE if TRMMdirName entered, 
							'CETRMMmin':floating-point representing the min rate in the CE if TRMMdirName entered}
	Assumptions::
	    Assumes we are dealing with MERG data which is 4kmx4km resolved, thus the smallest value 
        required according to Vila et al. (2008) is 2400km^2 
        therefore, 2400/16 = 150 contiguous squares
	'''
	
	frame = ma.empty((1,mergImgs.shape[1],mergImgs.shape[2]))
	CEcounter = 0
	frameCEcounter = 0
	frameNum = 0
	cloudElementEpsilon = 0.0
	cloudElementDict = {} 
	cloudElementCenter = []		#list with two elements [lat,lon] for the center of a CE
	prevFrameCEs = []			#list for CEs in previous frame
	currFrameCEs = []			#list for CEs in current frame
	cloudElementLat = []		#list for a particular CE's lat values
	cloudElementLon = []		#list for a particular CE's lon values
	cloudElementLatLons = []	#list for a particular CE's (lat,lon) values
	
	prevLatValue = 0.0
	prevLonValue = 0.0
	TIR_min = 0.0
	TIR_max = 0.0
	temporalRes = 3 # TRMM data is 3 hourly
	precipTotal = 0.0
	CETRMMList =[]
	precip =[]
	TRMMCloudElementLatLons =[]

	minCELatLimit = 0.0
	minCELonLimit = 0.0
	maxCELatLimit = 0.0
	maxCELonLimit = 0.0
	
	nygrd = len(LAT[:, 0]); nxgrd = len(LON[0, :])

	global MAINDIRECTORY
	MAINDIRECTORY = mainStrDir
	
	#openfile for storing ALL cloudElement information 
	cloudElementsFile = open((MAINDIRECTORY+'/textFiles/cloudElements.txt'),'wb')
	#openfile for storing cloudElement information meeting user criteria i.e. MCCs in this case
	cloudElementsUserFile = open((MAINDIRECTORY+'/textFiles/cloudElementsUserFile.txt'),'w')
	
	#NB in the TRMM files the info is hours since the time thus 00Z file has in 01, 02 and 03 times
	for t in xrange(mergImgs.shape[0]):
		#-------------------------------------------------
		# #textfile name for saving the data for arcgis
		# thisFileName = MAINDIRECTORY+'/' + (str(timelist[t])).replace(" ", "_") + '.txt'
		# cloudElementsTextFile = open(thisFileName,'w')
		#-------------------------------------------------

		#determine contiguous locations with temeperature below the warmest temp i.e. cloudElements in each frame
	   	frame, CEcounter = ndimage.measurements.label(mergImgs[t,:,:], structure=STRUCTURING_ELEMENT)
	   	frameCEcounter=0
		frameNum += 1

		#for each of the areas identified, check to determine if it a valid CE via an area and T requirement
	   	for count in xrange(CEcounter):
	   		#[0] is time dimension. Determine the actual values from the data
	   		#loc is a masked array
	   		try:
	   			loc = ndimage.find_objects(frame==(count+1))[0]
	   		except Exception, e:
	   			print "Error is ", e
	   			continue


	   		cloudElement = mergImgs[t,:,:][loc]
	   		labels, lcounter = ndimage.label(cloudElement)
	   		
	   		#determine the true lats and lons for this particular CE
   			cloudElementLat = LAT[loc[0],0]
   			cloudElementLon = LON[0,loc[1]] 
   			
	   		#determine number of boxes in this cloudelement
	   		numOfBoxes = np.count_nonzero(cloudElement)
	   		cloudElementArea = numOfBoxes*XRES*YRES

	   		#If the area is greater than the area required, or if the area is smaller than the suggested area, check if it meets a convective fraction requirement
	   		#consider as CE

	   		if cloudElementArea >= AREA_MIN or (cloudElementArea < AREA_MIN and ((ndimage.minimum(cloudElement, labels=labels))/float((ndimage.maximum(cloudElement, labels=labels)))) < CONVECTIVE_FRACTION ):

	   			#get some time information and labeling info
	   			frameTime = str(timelist[t])
	   			frameCEcounter +=1
	   			CEuniqueID = 'F'+str(frameNum)+'CE'+str(frameCEcounter) 

	   			#-------------------------------------------------
	    		#textfile name for accesing CE data using MATLAB code
				# thisFileName = MAINDIRECTORY+'/' + (str(timelist[t])).replace(" ", "_") + CEuniqueID +'.txt'
				# cloudElementsTextFile = open(thisFileName,'w')
				#-------------------------------------------------

				# ------ NETCDF File stuff for brightness temp stuff ------------------------------------
				thisFileName = MAINDIRECTORY +'/MERGnetcdfCEs/cloudElements'+ (str(timelist[t])).replace(" ", "_") + CEuniqueID +'.nc'
				currNetCDFCEData = Dataset(thisFileName, 'w', format='NETCDF4')
				currNetCDFCEData.description = 'Cloud Element '+CEuniqueID + ' temperature data'
				currNetCDFCEData.calendar = 'standard'
				currNetCDFCEData.conventions = 'COARDS'
				# dimensions
				currNetCDFCEData.createDimension('time', None)
				currNetCDFCEData.createDimension('lat', len(LAT[:,0]))
				currNetCDFCEData.createDimension('lon', len(LON[0,:]))
				# variables
				tempDims = ('time','lat', 'lon',)
				times = currNetCDFCEData.createVariable('time', 'f8', ('time',))
				times.units = 'hours since '+ str(timelist[t])[:-6]
				latitudes = currNetCDFCEData.createVariable('latitude', 'f8', ('lat',))
				longitudes = currNetCDFCEData.createVariable('longitude', 'f8', ('lon',))
				brightnesstemp = currNetCDFCEData.createVariable('brightnesstemp', 'i16',tempDims )
				brightnesstemp.units = 'Kelvin'
				# NETCDF data
				dates=[timelist[t]+timedelta(hours=0)]
				times[:] =  date2num(dates,units=times.units)
				longitudes[:] = LON[0,:]
				longitudes.units = "degrees_east" 
				longitudes.long_name = "Longitude" 

				latitudes[:] =  LAT[:,0]
				latitudes.units = "degrees_north"
				latitudes.long_name ="Latitude"
				
				#generate array of zeros for brightness temperature
				brightnesstemp1 = ma.zeros((1,len(latitudes), len(longitudes))).astype('int16')
				#-----------End most of NETCDF file stuff ------------------------------------

				#if other dataset (TRMM) assumed to be a precipitation dataset was entered
				if TRMMdirName:
					#------------------TRMM stuff -------------------------------------------------
					fileDate = ((str(timelist[t])).replace(" ", "")[:-8]).replace("-","")
					fileHr1 = (str(timelist[t])).replace(" ", "")[-8:-6]
					
					if int(fileHr1) % temporalRes == 0:
						fileHr = fileHr1
					else:
						fileHr = (int(fileHr1)/temporalRes) * temporalRes
					if fileHr < 10:
						fileHr = '0'+str(fileHr)
					else:
						str(fileHr)

					#open TRMM file for the resolution info and to create the appropriate sized grid
					TRMMfileName = TRMMdirName+'/3B42.'+ fileDate + "."+str(fileHr)+".7A.nc"
					
					TRMMData = Dataset(TRMMfileName,'r', format='NETCDF4')
					precipRate = TRMMData.variables['pcp'][:,:,:]
					latsrawTRMMData = TRMMData.variables['latitude'][:]
					lonsrawTRMMData = TRMMData.variables['longitude'][:]
					lonsrawTRMMData[lonsrawTRMMData > 180] = lonsrawTRMMData[lonsrawTRMMData>180] - 360.
					LONTRMM, LATTRMM = np.meshgrid(lonsrawTRMMData, latsrawTRMMData)

					nygrdTRMM = len(LATTRMM[:,0]); nxgrdTRMM = len(LONTRMM[0,:])
					precipRateMasked = ma.masked_array(precipRate, mask=(precipRate < 0.0))
					#---------regrid the TRMM data to the MERG dataset ----------------------------------
					#regrid using the do_regrid stuff from the Apache OCW 
					regriddedTRMM = ma.zeros((0, nygrd, nxgrd))
					regriddedTRMM = utils.doRegrid(precipRateMasked[0,:,:], LATTRMM,  LONTRMM, LAT, LON, order=1, mdi= -999999999)
					#----------------------------------------------------------------------------------
		
					# #get the lat/lon info from cloudElement
					#get the lat/lon info from the file
					latCEStart = LAT[0][0]
					latCEEnd = LAT[-1][0]
					lonCEStart = LON[0][0]
					lonCEEnd = LON[0][-1]
					
					#get the lat/lon info for TRMM data (different resolution)
					latStartT = utils.findNearest(latsrawTRMMData, latCEStart)
					latEndT = utils.findNearest(latsrawTRMMData, latCEEnd)
					lonStartT = utils.findNearest(lonsrawTRMMData, lonCEStart)
					lonEndT = utils.findNearest(lonsrawTRMMData, lonCEEnd)
					latStartIndex = np.where(latsrawTRMMData == latStartT)
					latEndIndex = np.where(latsrawTRMMData == latEndT)
					lonStartIndex = np.where(lonsrawTRMMData == lonStartT)
					lonEndIndex = np.where(lonsrawTRMMData == lonEndT)

					#get the relevant TRMM info 
					CEprecipRate = precipRate[:,(latStartIndex[0][0]-1):latEndIndex[0][0],(lonStartIndex[0][0]-1):lonEndIndex[0][0]]
					TRMMData.close()
					
					# ------ NETCDF File info for writing TRMM CE rainfall ------------------------------------
					thisFileName = MAINDIRECTORY+'/TRMMnetcdfCEs/TRMM' + (str(timelist[t])).replace(" ", "_") + CEuniqueID +'.nc'
					currNetCDFTRMMData = Dataset(thisFileName, 'w', format='NETCDF4')
					currNetCDFTRMMData.description = 'Cloud Element '+CEuniqueID + ' precipitation data'
					currNetCDFTRMMData.calendar = 'standard'
					currNetCDFTRMMData.conventions = 'COARDS'
					# dimensions
					currNetCDFTRMMData.createDimension('time', None)
					currNetCDFTRMMData.createDimension('lat', len(LAT[:,0]))
					currNetCDFTRMMData.createDimension('lon', len(LON[0,:]))
					
					# variables
					TRMMprecip = ('time','lat', 'lon',)
					times = currNetCDFTRMMData.createVariable('time', 'f8', ('time',))
					times.units = 'hours since '+ str(timelist[t])[:-6]
					latitude = currNetCDFTRMMData.createVariable('latitude', 'f8', ('lat',))
					longitude = currNetCDFTRMMData.createVariable('longitude', 'f8', ('lon',))
					rainFallacc = currNetCDFTRMMData.createVariable('precipitation_Accumulation', 'f8',TRMMprecip )
					rainFallacc.units = 'mm'

					longitude[:] = LON[0,:]
					longitude.units = "degrees_east" 
					longitude.long_name = "Longitude" 

					latitude[:] =  LAT[:,0]
					latitude.units = "degrees_north"
					latitude.long_name ="Latitude"

					finalCETRMMvalues = ma.zeros((brightnesstemp.shape))
					#-----------End most of NETCDF file stuff ------------------------------------

	   			#populate cloudElementLatLons by unpacking the original values from loc to get the actual value for lat and lon
    			#TODO: KDW - too dirty... play with itertools.izip or zip and the enumerate with this
    			# 			as cloudElement is masked
				for index,value in np.ndenumerate(cloudElement):
					if value != 0 : 
						lat_index,lon_index = index
						lat_lon_tuple = (cloudElementLat[lat_index], cloudElementLon[lon_index],value)

						#generate the comma separated file for GIS
						cloudElementLatLons.append(lat_lon_tuple)

						#temp data for CE NETCDF file
						brightnesstemp1[0,int(np.where(LAT[:,0]==cloudElementLat[lat_index])[0]),int(np.where(LON[0,:]==cloudElementLon[lon_index])[0])] = value
						
						if TRMMdirName:
							finalCETRMMvalues[0,int(np.where(LAT[:,0]==cloudElementLat[lat_index])[0]),int(np.where(LON[0,:]==cloudElementLon[lon_index])[0])] = regriddedTRMM[int(np.where(LAT[:,0]==cloudElementLat[lat_index])[0]),int(np.where(LON[0,:]==cloudElementLon[lon_index])[0])]
							CETRMMList.append((cloudElementLat[lat_index], cloudElementLon[lon_index], finalCETRMMvalues[0,cloudElementLat[lat_index], cloudElementLon[lon_index]]))


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
					TRMMArea = TRMMnumOfBoxes*XRES*YRES
					try:
						maxCEprecipRate = np.max(finalCETRMMvalues[np.nonzero(finalCETRMMvalues)])
						minCEprecipRate = np.min(finalCETRMMvalues[np.nonzero(finalCETRMMvalues)])
					except:
						pass

				#sort cloudElementLatLons by lats
				cloudElementLatLons.sort(key=lambda tup: tup[0])	

				#determine if the cloud element the shape 
				cloudElementEpsilon = eccentricity (cloudElement)
	   			cloudElementsUserFile.write("\n\nTime is: %s" %(str(timelist[t])))
	   			cloudElementsUserFile.write("\nCEuniqueID is: %s" %CEuniqueID)
	   			latCenter, lonCenter = ndimage.measurements.center_of_mass(cloudElement, labels=labels)
	   			
	   			#latCenter and lonCenter are given according to the particular array defining this CE
	   			#so you need to convert this value to the overall domain truth
	   			latCenter = cloudElementLat[round(latCenter)]
	   			lonCenter = cloudElementLon[round(lonCenter)]
	   			cloudElementsUserFile.write("\nCenter (lat,lon) is: %.2f\t%.2f" %(latCenter, lonCenter))
	   			cloudElementCenter.append(latCenter)
	   			cloudElementCenter.append(lonCenter)
	   			cloudElementsUserFile.write("\nNumber of boxes are: %d" %numOfBoxes)
	   			cloudElementsUserFile.write("\nArea is: %.4f km^2" %(cloudElementArea))
				cloudElementsUserFile.write("\nAverage brightness temperature is: %.4f K" %ndimage.mean(cloudElement, labels=labels))
				cloudElementsUserFile.write("\nMin brightness temperature is: %.4f K" %ndimage.minimum(cloudElement, labels=labels))
				cloudElementsUserFile.write("\nMax brightness temperature is: %.4f K" %ndimage.maximum(cloudElement, labels=labels))
				cloudElementsUserFile.write("\nBrightness temperature variance is: %.4f K" %ndimage.variance(cloudElement, labels=labels))
				cloudElementsUserFile.write("\nConvective fraction is: %.4f " %(((ndimage.minimum(cloudElement, labels=labels))/float((ndimage.maximum(cloudElement, labels=labels))))*100.0))
				cloudElementsUserFile.write("\nEccentricity is: %.4f " %(cloudElementEpsilon))
				#populate the dictionary
				if TRMMdirName:
					cloudElementDict = {'uniqueID': CEuniqueID, 'cloudElementTime': timelist[t],'cloudElementLatLon': cloudElementLatLons, 'cloudElementCenter':cloudElementCenter, 'cloudElementArea':cloudElementArea, 'cloudElementEccentricity':cloudElementEpsilon, 'cloudElementTmax':TIR_max, 'cloudElementTmin': TIR_min, 'cloudElementPrecipTotal':precipTotal,'cloudElementLatLonTRMM':CETRMMList, 'TRMMArea': TRMMArea,'CETRMMmax':maxCEprecipRate, 'CETRMMmin':minCEprecipRate}
				else:
					cloudElementDict = {'uniqueID': CEuniqueID, 'cloudElementTime': timelist[t],'cloudElementLatLon': cloudElementLatLons, 'cloudElementCenter':cloudElementCenter, 'cloudElementArea':cloudElementArea, 'cloudElementEccentricity':cloudElementEpsilon, 'cloudElementTmax':TIR_max, 'cloudElementTmin': TIR_min,}
				
				#current frame list of CEs
				currFrameCEs.append(cloudElementDict)
				
				#draw the graph node
				CLOUD_ELEMENT_GRAPH.add_node(CEuniqueID, cloudElementDict)
				
				if frameNum != 1:
					for cloudElementDict in prevFrameCEs:
						thisCElen = len(cloudElementLatLons)
						percentageOverlap, areaOverlap = cloudElementOverlap(cloudElementLatLons, cloudElementDict['cloudElementLatLon'])
						
						#change weights to integers because the built in shortest path chokes on floating pts according to Networkx doc
						#according to Goyens et al, two CEs are considered related if there is atleast 95% overlap between them for consecutive imgs a max of 2 hrs apart
						if percentageOverlap >= 0.95: 
							CLOUD_ELEMENT_GRAPH.add_edge(cloudElementDict['uniqueID'], CEuniqueID, weight=edgeWeight[0])
							
						elif percentageOverlap >= 0.90 and percentageOverlap < 0.95 :
							CLOUD_ELEMENT_GRAPH.add_edge(cloudElementDict['uniqueID'], CEuniqueID, weight=edgeWeight[1])

						elif areaOverlap >= MIN_OVERLAP:
							CLOUD_ELEMENT_GRAPH.add_edge(cloudElementDict['uniqueID'], CEuniqueID, weight=edgeWeight[2])

    			else:
    				#TODO: remove this else as we only wish for the CE details
    				#ensure only the non-zero elements are considered
    				#store intel in allCE file
    				labels, _ = ndimage.label(cloudElement)
    				cloudElementsFile.write("\n-----------------------------------------------")
    				cloudElementsFile.write("\n\nTime is: %s" %(str(timelist[t])))
    				# cloudElementLat = LAT[loc[0],0]
    				# cloudElementLon = LON[0,loc[1]] 
    				
    				#populate cloudElementLatLons by unpacking the original values from loc
    				#TODO: KDW - too dirty... play with itertools.izip or zip and the enumerate with this
    				# 			as cloudElement is masked
    				for index,value in np.ndenumerate(cloudElement):
    					if value != 0 : 
    						lat_index,lon_index = index
    						lat_lon_tuple = (cloudElementLat[lat_index], cloudElementLon[lon_index])
    						cloudElementLatLons.append(lat_lon_tuple)
	
    				cloudElementsFile.write("\nLocation of rejected CE (lat,lon) points are: %s" %cloudElementLatLons)
    				#latCenter and lonCenter are given according to the particular array defining this CE
		   			#so you need to convert this value to the overall domain truth
    				latCenter, lonCenter = ndimage.measurements.center_of_mass(cloudElement, labels=labels)
    				latCenter = cloudElementLat[round(latCenter)]
    				lonCenter = cloudElementLon[round(lonCenter)]
    				cloudElementsFile.write("\nCenter (lat,lon) is: %.2f\t%.2f" %(latCenter, lonCenter))
    				cloudElementsFile.write("\nNumber of boxes are: %d" %numOfBoxes)
    				cloudElementsFile.write("\nArea is: %.4f km^2" %(cloudElementArea))
    				cloudElementsFile.write("\nAverage brightness temperature is: %.4f K" %ndimage.mean(cloudElement, labels=labels))
    				cloudElementsFile.write("\nMin brightness temperature is: %.4f K" %ndimage.minimum(cloudElement, labels=labels))
    				cloudElementsFile.write("\nMax brightness temperature is: %.4f K" %ndimage.maximum(cloudElement, labels=labels))
    				cloudElementsFile.write("\nBrightness temperature variance is: %.4f K" %ndimage.variance(cloudElement, labels=labels))
    				cloudElementsFile.write("\nConvective fraction is: %.4f " %(((ndimage.minimum(cloudElement, labels=labels))/float((ndimage.maximum(cloudElement, labels=labels))))*100.0))
    				cloudElementsFile.write("\nEccentricity is: %.4f " %(cloudElementEpsilon))
    				cloudElementsFile.write("\n-----------------------------------------------")
    				
			#reset list for the next CE
			nodeExist = False
			cloudElementCenter=[]
			cloudElement = []
			cloudElementLat=[]
			cloudElementLon =[]
			cloudElementLatLons =[]
			brightnesstemp1 =[]
			brightnesstemp =[]
			finalCETRMMvalues =[]
			CEprecipRate =[]
			CETRMMList =[]
			precipTotal = 0.0
			precip=[]
			TRMMCloudElementLatLons=[]
			
		#reset for the next time
		prevFrameCEs =[]
		prevFrameCEs = currFrameCEs
		currFrameCEs =[]
		    			
	cloudElementsFile.close
	cloudElementsUserFile.close
	#if using ARCGIS data store code, uncomment this file close line
	#cloudElementsTextFile.close

	#clean up graph - remove parent and childless nodes
	outAndInDeg = CLOUD_ELEMENT_GRAPH.degree_iter()
	toRemove = [node[0] for node in outAndInDeg if node[1]<1]
	CLOUD_ELEMENT_GRAPH.remove_nodes_from(toRemove)
	
	print "number of nodes are: ", CLOUD_ELEMENT_GRAPH.number_of_nodes()
	print "number of edges are: ", CLOUD_ELEMENT_GRAPH.number_of_edges()
	print ("*"*80)

	#hierachial graph output
	graphTitle = "Cloud Elements observed over somewhere from 0000Z to 0000Z" 
	plotting.drawGraph(CLOUD_ELEMENT_GRAPH, graphTitle, MAINDIRECTORY, edgeWeight)

	return CLOUD_ELEMENT_GRAPH	
#******************************************************************
def findPrecipRate(TRMMdirName, timelist):
	''' 
	Purpose:: 
		Determines the precipitation rates for MCSs found if TRMMdirName was not entered in findCloudElements this can be used

	Input:: 
		TRMMdirName: a string representing the directory for the original TRMM netCDF files
		timelist: a list of python datatimes

    Output:: a list of dictionary of the TRMM data 
    	NB: also creates netCDF with TRMM data for each CE (for post processing) index
    		in MAINDIRECTORY/TRMMnetcdfCEs
   
    Assumptions:: Assumes that findCloudElements was run without the TRMMdirName value 
 
	'''
	allCEnodesTRMMdata =[]
	TRMMdataDict={}
	precipTotal = 0.0

	os.chdir((MAINDIRECTORY+'/MERGnetcdfCEs/'))
	imgFilename = ''
	temporalRes = 3 #3 hours for TRMM
	
	#sort files
	files = filter(os.path.isfile, glob.glob("*.nc"))
	files.sort(key=lambda x: os.path.getmtime(x))
	
	for afile in files:
		fullFname = os.path.splitext(afile)[0]
		noFrameExtension = (fullFname.replace("_","")).split('F')[0]
		CEuniqueID = 'F' +(fullFname.replace("_","")).split('F')[1]
		fileDateTimeChar = (noFrameExtension.replace(":","")).split('s')[1]
		fileDateTime = fileDateTimeChar.replace("-","")
		fileDate = fileDateTime[:-6]
		fileHr1=fileDateTime[-6:-4]

		cloudElementData = Dataset(afile,'r', format='NETCDF4')
		brightnesstemp1 = cloudElementData.variables['brightnesstemp'][:,:,:] 
		latsrawCloudElements = cloudElementData.variables['latitude'][:]
		lonsrawCloudElements = cloudElementData.variables['longitude'][:]
		
		brightnesstemp = np.squeeze(brightnesstemp1, axis=0)
		
		if int(fileHr1) % temporalRes == 0:
			fileHr = fileHr1
		else:
			fileHr = (int(fileHr1)/temporalRes) * temporalRes
		
		if fileHr < 10:
			fileHr = '0'+str(fileHr)
		else:
			str(fileHr)

		TRMMfileName = TRMMdirName+"/3B42."+ str(fileDate) + "."+str(fileHr)+".7A.nc"
		TRMMData = Dataset(TRMMfileName,'r', format='NETCDF4')
		precipRate = TRMMData.variables['pcp'][:,:,:]
		latsrawTRMMData = TRMMData.variables['latitude'][:]
		lonsrawTRMMData = TRMMData.variables['longitude'][:]
		lonsrawTRMMData[lonsrawTRMMData > 180] = lonsrawTRMMData[lonsrawTRMMData>180] - 360.
		LONTRMM, LATTRMM = np.meshgrid(lonsrawTRMMData, latsrawTRMMData)

		#nygrdTRMM = len(LATTRMM[:,0]); nxgrd = len(LONTRMM[0,:])
		nygrd = len(LAT[:, 0]); nxgrd = len(LON[0, :])

		precipRateMasked = ma.masked_array(precipRate, mask=(precipRate < 0.0))
		#---------regrid the TRMM data to the MERG dataset ----------------------------------
		#regrid using the do_regrid stuff from the Apache OCW 
		regriddedTRMM = ma.zeros((0, nygrd, nxgrd))
		regriddedTRMM = utils.doRegrid(precipRateMasked[0,:,:], LATTRMM,  LONTRMM, LAT, LON, order=1, mdi= -999999999)
		#----------------------------------------------------------------------------------

		# #get the lat/lon info from
		latCEStart = LAT[0][0]
		latCEEnd = LAT[-1][0]
		lonCEStart = LON[0][0]
		lonCEEnd = LON[0][-1]

		#get the lat/lon info for TRMM data (different resolution)
		latStartT = utils.findNearest(latsrawTRMMData, latCEStart)
		latEndT = utils.findNearest(latsrawTRMMData, latCEEnd)
		lonStartT = utils.findNearest(lonsrawTRMMData, lonCEStart)
		lonEndT = utils.findNearest(lonsrawTRMMData, lonCEEnd)
		latStartIndex = np.where(latsrawTRMMData == latStartT)
		latEndIndex = np.where(latsrawTRMMData == latEndT)
		lonStartIndex = np.where(lonsrawTRMMData == lonStartT)
		lonEndIndex = np.where(lonsrawTRMMData == lonEndT)

		#get the relevant TRMM info 
		CEprecipRate = precipRate[:,(latStartIndex[0][0]-1):latEndIndex[0][0],(lonStartIndex[0][0]-1):lonEndIndex[0][0]]
		TRMMData.close()
			
		
		# ------ NETCDF File stuff ------------------------------------
		thisFileName = MAINDIRECTORY+'/TRMMnetcdfCEs/'+ fileDateTime + CEuniqueID +'.nc'
		currNetCDFTRMMData = Dataset(thisFileName, 'w', format='NETCDF4')
		currNetCDFTRMMData.description = 'Cloud Element '+CEuniqueID + ' rainfall data'
		currNetCDFTRMMData.calendar = 'standard'
		currNetCDFTRMMData.conventions = 'COARDS'
		# dimensions
		currNetCDFTRMMData.createDimension('time', None)
		currNetCDFTRMMData.createDimension('lat', len(LAT[:,0]))
		currNetCDFTRMMData.createDimension('lon', len(LON[0,:]))
		# variables
		TRMMprecip = ('time','lat', 'lon',)
		times = currNetCDFTRMMData.createVariable('time', 'f8', ('time',))
		times.units = 'hours since '+ fileDateTime[:-6] 
		latitude = currNetCDFTRMMData.createVariable('latitude', 'f8', ('lat',))
		longitude = currNetCDFTRMMData.createVariable('longitude', 'f8', ('lon',))
		rainFallacc = currNetCDFTRMMData.createVariable('precipitation_Accumulation', 'f8',TRMMprecip )
		rainFallacc.units = 'mm'

		longitude[:] = LON[0,:]
		longitude.units = "degrees_east" 
		longitude.long_name = "Longitude" 

		latitude[:] =  LAT[:,0]
		latitude.units = "degrees_north"
		latitude.long_name ="Latitude"

		finalCETRMMvalues = ma.zeros((brightnesstemp1.shape))
		#-----------End most of NETCDF file stuff ------------------------------------	
		for index,value in np.ndenumerate(brightnesstemp):
			lat_index, lon_index = index
			currTimeValue = 0
			if value > 0:

				finalCETRMMvalues[0,lat_index,lon_index] = regriddedTRMM[int(np.where(LAT[:,0]==LAT[lat_index,0])[0]), int(np.where(LON[0,:]==LON[0,lon_index])[0])]
				

		rainFallacc[:] = finalCETRMMvalues
		currNetCDFTRMMData.close()

		for index, value in np.ndenumerate(finalCETRMMvalues):
			precipTotal += value 

		TRMMnumOfBoxes = np.count_nonzero(finalCETRMMvalues)
		TRMMArea = TRMMnumOfBoxes*XRES*YRES	

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
		for eachdict in CLOUD_ELEMENT_GRAPH.nodes(CEuniqueID):
			if eachdict[1]['uniqueID'] == CEuniqueID:
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
		latsrawTRMMData =[]
		lonsrawTRMMData = []
		latsrawCloudElements=[]
		lonsrawCloudElements=[]
		finalCETRMMvalues =[]
		CEprecipRate =[]
		brightnesstemp =[]
		TRMMdataDict ={}

	return allCEnodesTRMMdata
#******************************************************************	
def findCloudClusters(CEGraph):
	'''
	Purpose:: 
		Determines the cloud clusters properties from the subgraphs in 
	    the graph i.e. prunes the graph according to the minimum depth

	Input:: 
		CEGraph: a Networkx directed graph of the CEs with weighted edges
		according the area overlap between nodes (CEs) of consectuive frames
    
    Output:: 
    	PRUNED_GRAPH: a Networkx directed graph of with CCs/ MCSs

	'''

	seenNode = []
	allMCSLists =[]
	pathDictList =[]
	pathList=[]

	cloudClustersFile = open((MAINDIRECTORY+'/textFiles/cloudClusters.txt'),'wb')
	
	for eachNode in CEGraph:
		#check if the node has been seen before
		if eachNode not in dict(enumerate(zip(*seenNode))):
			#look for all trees associated with node as the root
			thisPathDistanceAndLength = nx.single_source_dijkstra(CEGraph, eachNode)
			#determine the actual shortestPath and minimum depth/length
			maxDepthAndMinPath = findMaxDepthAndMinPath(thisPathDistanceAndLength)
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
				cloudClustersFile.write("\nSubtree pathlength is %d and path is %s" %(maxPathLength, shortestPath))
				#update seenNode info
				seenNode.append(shortestPath)	

	print "pruned graph"
	print "number of nodes are: ", PRUNED_GRAPH.number_of_nodes()
	print "number of edges are: ", PRUNED_GRAPH.number_of_edges()
	print ("*"*80)		
					
	graphTitle = "Cloud Clusters observed over somewhere during sometime"
	plotting.drawGraph(PRUNED_GRAPH, graphTitle, MAINDIRECTORY, edgeWeight)
	cloudClustersFile.close
	
	return PRUNED_GRAPH  
#******************************************************************
def findMCC (prunedGraph):
	'''
	Purpose:: 
		Determines if subtree is a MCC according to Laurent et al 1998 criteria

	Input:: 
		prunedGraph: a Networkx Graph representing the CCs 

    Output:: 
    	finalMCCList: a list of list of tuples representing a MCC
             
    Assumptions: 
    	frames are ordered and are equally distributed in time e.g. hrly satellite images
 
	'''
	MCCList = []
	MCSList = []
	definiteMCC = []
	definiteMCS = []
	eachList =[]
	eachMCCList =[]
	maturing = False
	decaying = False
	fNode = ''
	lNode = ''
	removeList =[]
	imgCount = 0
	imgTitle =''
	
	maxShieldNode = ''
	orderedPath =[]
	treeTraversalList =[]
	definiteMCCFlag = False
	unDirGraph = nx.Graph()
	aSubGraph = nx.DiGraph()
	definiteMCSFlag = False

	
	#connected_components is not available for DiGraph, so generate graph as undirected 
	unDirGraph = PRUNED_GRAPH.to_undirected()
	subGraph = nx.connected_component_subgraphs(unDirGraph)

	#for each path in the subgraphs determined
	for path in subGraph:
		#definite is a subTree provided the duration is longer than 3 hours

		if len(path.nodes()) > MIN_MCS_DURATION:
			orderedPath = path.nodes()
			orderedPath.sort(key=lambda item:(len(item.split('C')[0]), item.split('C')[0]))
			#definiteMCS.append(orderedPath)

			#build back DiGraph for checking purposes/paper purposes
			aSubGraph.add_nodes_from(path.nodes())	
			for eachNode in path.nodes():
				if prunedGraph.predecessors(eachNode):
					for node in prunedGraph.predecessors(eachNode):
						aSubGraph.add_edge(node,eachNode,weight=edgeWeight[0])

				if prunedGraph.successors(eachNode):
					for node in prunedGraph.successors(eachNode):
						aSubGraph.add_edge(eachNode,node,weight=edgeWeight[0])
			imgTitle = 'CC'+str(imgCount+1)
			plotting.drawGraph(aSubGraph, imgTitle, MAINDIRECTORY, edgeWeight) #for eachNode in path:
			imgCount +=1
			#----------end build back ---------------------------------------------

			mergeList, splitList = hasMergesOrSplits(path)	
			#add node behavior regarding neutral, merge, split or both
			for node in path:
				if node in mergeList and node in splitList:
					addNodeBehaviorIdentifier(node,'B')
				elif node in mergeList and not node in splitList:
					addNodeBehaviorIdentifier(node,'M')
				elif node in splitList and not node in mergeList:
					addNodeBehaviorIdentifier(node,'S')
				else:
					addNodeBehaviorIdentifier(node,'N')
			

			#Do the first part of checking for the MCC feature
			#find the path
			treeTraversalList = traverseTree(aSubGraph, orderedPath[0],[],[])
			#print "treeTraversalList is ", treeTraversalList
			#check the nodes to determine if a MCC on just the area criteria (consecutive nodes meeting the area and temp requirements)
			MCCList = checkedNodesMCC(prunedGraph, treeTraversalList)
			for aDict in MCCList:
				for eachNode in aDict["fullMCSMCC"]:
					addNodeMCSIdentifier(eachNode[0],eachNode[1])
				
			#do check for if MCCs overlap
			if MCCList:
				if len(MCCList) > 1:
					for count in range(len(MCCList)): #for eachDict in MCCList:
						#if there are more than two lists
						if count >= 1:
							#and the first node in this list
							eachList = list(x[0] for x in MCCList[count]["possMCCList"])
							eachList.sort(key=lambda nodeID:(len(nodeID.split('C')[0]), nodeID.split('C')[0]))
							if eachList:
								fNode = eachList[0]
								#get the lastNode in the previous possMCC list
								eachList = list(x[0] for x in MCCList[(count-1)]["possMCCList"])
								eachList.sort(key=lambda nodeID:(len(nodeID.split('C')[0]), nodeID.split('C')[0]))
								if eachList:
									lNode = eachList[-1]
									if lNode in CLOUD_ELEMENT_GRAPH.predecessors(fNode):
										for aNode in CLOUD_ELEMENT_GRAPH.predecessors(fNode):
											if aNode in eachList and aNode == lNode:
												#if edge_data is equal or less than to the exisitng edge in the tree append one to the other
												if CLOUD_ELEMENT_GRAPH.get_edge_data(aNode,fNode)['weight'] <= CLOUD_ELEMENT_GRAPH.get_edge_data(lNode,fNode)['weight']:
													MCCList[count-1]["possMCCList"].extend(MCCList[count]["possMCCList"]) 
													MCCList[count-1]["fullMCSMCC"].extend(MCCList[count]["fullMCSMCC"])
													MCCList[count-1]["durationAandB"] +=  MCCList[count]["durationAandB"]
													MCCList[count-1]["CounterCriteriaA"] += MCCList[count]["CounterCriteriaA"]
													MCCList[count-1]["highestMCCnode"] = MCCList[count]["highestMCCnode"]
													MCCList[count-1]["frameNum"] = MCCList[count]["frameNum"] 
													removeList.append(count)
				#update the MCCList
				if removeList:
					for i in removeList:
						if (len(MCCList)-1) > i:
							del MCCList[i]
							removeList =[]
				
			#check if the nodes also meet the duration criteria and the shape crieria
			for eachDict in MCCList:
				#order the fullMCSMCC list, then run maximum extent and eccentricity criteria 
				if (eachDict["durationAandB"] * TRES) >= MINIMUM_DURATION and (eachDict["durationAandB"] * TRES) <= MAXIMUM_DURATION:
					eachList = list(x[0] for x in eachDict["fullMCSMCC"])
					eachList.sort(key=lambda nodeID:(len(nodeID.split('C')[0]), nodeID.split('C')[0]))
					eachMCCList = list(x[0] for x in eachDict["possMCCList"])
					eachMCCList.sort(key=lambda nodeID:(len(nodeID.split('C')[0]), nodeID.split('C')[0]))
					
					#update the nodemcsidentifer behavior
					#find the first element eachMCCList in eachList, and ensure everything ahead of it is indicated as 'I', 
					#find last element in eachMCCList in eachList and ensure everything after it is indicated as 'D'
					#ensure that everything between is listed as 'M'
					for eachNode in eachList[:(eachList.index(eachMCCList[0]))]: 
						addNodeMCSIdentifier(eachNode,'I')

					addNodeMCSIdentifier(eachMCCList[0],'M')

					for eachNode in eachList[(eachList.index(eachMCCList[-1])+1):]:
						addNodeMCSIdentifier(eachNode, 'D')

					#update definiteMCS list
					for eachNode in orderedPath[(orderedPath.index(eachMCCList[-1])+1):]:
						addNodeMCSIdentifier(eachNode, 'D')

					#run maximum extent and eccentricity criteria
					maxExtentNode, definiteMCCFlag = maxExtentAndEccentricity(eachList)
					#print "maxExtentNode, definiteMCCFlag ", maxExtentNode, definiteMCCFlag
					if definiteMCCFlag == True:
						definiteMCC.append(eachList)


			definiteMCS.append(orderedPath)
			
			#reset for next subGraph	
			aSubGraph.clear()
			orderedPath=[]
			MCCList =[]
			MCSList =[]
			definiteMCSFlag = False
		
	return definiteMCC, definiteMCS
#******************************************************************
def traverseTree(subGraph,node, stack, checkedNodes=None):
	'''
	Purpose:: 
		To traverse a tree using a modified depth-first iterative deepening (DFID) search algorithm 

	Input:: 
		subGraph: a Networkx DiGraph representing a CC
			lengthOfsubGraph: an integer representing the length of the subgraph
			node: a string representing the node currently being checked
			stack: a list of strings representing a list of nodes in a stack functionality 
					i.e. Last-In-First-Out (LIFO) for sorting the information from each visited node
			checkedNodes: a list of strings representing the list of the nodes in the traversal
    
    Output:: 
    	checkedNodes: a list of strings representing the list of the nodes in the traversal

    Assumptions: 
    	frames are ordered and are equally distributed in time e.g. hrly satellite images
 
	'''
	if len(checkedNodes) == len(subGraph):
		return checkedNodes

	if not checkedNodes:
		stack =[]
		checkedNodes.append(node)
		
	#check one level infront first...if something does exisit, stick it at the front of the stack
	upOneLevel = subGraph.predecessors(node)
	downOneLevel = subGraph.successors(node)
	for parent in upOneLevel:
		if parent not in checkedNodes and parent not in stack:
			for child in downOneLevel:
				if child not in checkedNodes and child not in stack:
					stack.insert(0,child)
		
			stack.insert(0,parent)	

	for child in downOneLevel:
		if child not in checkedNodes and child not in stack:
			if len(subGraph.predecessors(child)) > 1 or node in checkedNodes:
				stack.insert(0,child)
			else:
				stack.append(child)		
	
	for eachNode in stack:
		if eachNode not in checkedNodes:
			checkedNodes.append(eachNode)
			return traverseTree(subGraph, eachNode, stack, checkedNodes)
	
	return checkedNodes 
#******************************************************************
def checkedNodesMCC (prunedGraph, nodeList):
	'''
	Purpose :: 
		Determine if this path is (or is part of) a MCC and provides 
	    preliminary information regarding the stages of the feature

	Input:: 
		prunedGraph: a Networkx Graph representing all the cloud clusters 
		nodeList: list of strings (CE ID) from the traversal
		
	Output:: 
		potentialMCCList: list of dictionaries representing all possible MCC within the path
			dictionary = {"possMCCList":[(node,'I')], "fullMCSMCC":[(node,'I')], "CounterCriteriaA": CounterCriteriaA, "durationAandB": durationAandB}
	'''
	
	CounterCriteriaAFlag = False
	CounterCriteriaBFlag = False
	INITIATIONFLAG = False
	MATURITYFLAG = False
	DECAYFLAG = False
	thisdict = {} #will have the same items as the cloudElementDict 
	cloudElementAreaB = 0.0
	cloudElementAreaA = 0.0
	epsilon = 0.0
	frameNum =0
	oldNode =''
	potentialMCCList =[]
	durationAandB = 0

	#check for if the list contains only one string/node
	if type(nodeList) is str:
		oldNode=nodeList
		nodeList =[]
		nodeList.append(oldNode)

	for node in nodeList:
		thisdict = thisDict(node)
		CounterCriteriaAFlag = False
		CounterCriteriaBFlag = False
		existingFrameFlag = False

		if thisdict['cloudElementArea'] >= OUTER_CLOUD_SHIELD_AREA:
			CounterCriteriaAFlag = True
			INITIATIONFLAG = True
			MATURITYFLAG = False

			#check if criteriaA is met
			cloudElementAreaA, criteriaA = checkCriteria(thisdict['cloudElementLatLon'], OUTER_CLOUD_SHIELD_TEMPERATURE)
			#TODO: calcuate the eccentricity at this point and read over????or create a new field in the dict
			
			if cloudElementAreaA >= OUTER_CLOUD_SHIELD_AREA:
				#check if criteriaB is met
				cloudElementAreaB,criteriaB = checkCriteria(thisdict['cloudElementLatLon'], INNER_CLOUD_SHIELD_TEMPERATURE)
				
				#if Criteria A and B have been met, then the MCC is initiated, i.e. store node as potentialMCC
		   		if cloudElementAreaB >= INNER_CLOUD_SHIELD_AREA:
		   			#TODO: add another field to the dictionary for the OUTER_AREA_SHIELD area
		   			CounterCriteriaBFlag = True
		   			#append this information on to the dictionary
		   			addInfothisDict(node, cloudElementAreaB, criteriaB)
		   			INITIATIONFLAG = False
		   			MATURITYFLAG = True
		   			stage = 'M'
		   			potentialMCCList = updateMCCList(prunedGraph, potentialMCCList,node,stage, CounterCriteriaAFlag, CounterCriteriaBFlag) 			
	   			else:
	   				#criteria B failed
	   				CounterCriteriaBFlag = False
	   				if INITIATIONFLAG == True:
	   					stage = 'I'   					
	   					potentialMCCList = updateMCCList(prunedGraph, potentialMCCList,node,stage, CounterCriteriaAFlag, CounterCriteriaBFlag)

	   				elif (INITIATIONFLAG == False and MATURITYFLAG == True) or DECAYFLAG==True:
	   					DECAYFLAG = True
	   					MATURITYFLAG = False
	   					stage = 'D'
	   					potentialMCCList = updateMCCList(prunedGraph, potentialMCCList,node,stage, CounterCriteriaAFlag, CounterCriteriaBFlag)
	   		else:
	   			#criteria A failed
	   			CounterCriteriaAFlag = False
	   			CounterCriteriaBFlag = False
	   			#add as a CE before or after the main feature
				if INITIATIONFLAG == True or (INITIATIONFLAG == False and MATURITYFLAG == True):
					stage ="I"
					potentialMCCList = updateMCCList(prunedGraph, potentialMCCList,node,stage, CounterCriteriaAFlag, CounterCriteriaBFlag)
	   			elif (INITIATIONFLAG == False and MATURITYFLAG == False) or DECAYFLAG == True:
	   				stage = "D"
	   				DECAYFLAG = True
	   				potentialMCCList = updateMCCList(prunedGraph, potentialMCCList,node,stage, CounterCriteriaAFlag, CounterCriteriaBFlag)
	   			elif (INITIATIONFLAG == False and MATURITYFLAG == False and DECAYFLAG == False):
	   				stage ="I"
	   				potentialMCCList = updateMCCList(prunedGraph, potentialMCCList,node,stage, CounterCriteriaAFlag, CounterCriteriaBFlag)


   		else:
   			#criteria A failed
   			CounterCriteriaAFlag = False
   			CounterCriteriaBFlag = False
   			#add as a CE before or after the main feature
			if INITIATIONFLAG == True or (INITIATIONFLAG == False and MATURITYFLAG == True):
				stage ="I"
				potentialMCCList = updateMCCList(prunedGraph, potentialMCCList,node,stage, CounterCriteriaAFlag, CounterCriteriaBFlag)
   			elif (INITIATIONFLAG == False and MATURITYFLAG == False) or DECAYFLAG == True:
   				stage = "D"
   				DECAYFLAG = True
   				potentialMCCList = updateMCCList(prunedGraph, potentialMCCList,node,stage, CounterCriteriaAFlag, CounterCriteriaBFlag)
   			elif (INITIATIONFLAG == False and MATURITYFLAG == False and DECAYFLAG == False):
   				stage ="I"
   				potentialMCCList = updateMCCList(prunedGraph, potentialMCCList,node,stage, CounterCriteriaAFlag, CounterCriteriaBFlag)

	return potentialMCCList
#******************************************************************
def updateMCCList(prunedGraph, potentialMCCList,node,stage, CounterCriteriaAFlag, CounterCriteriaBFlag):
	'''
	Purpose:: 
		Utility function to determine if a path is (or is part of) a MCC and provides 
	           preliminary information regarding the stages of the feature

	Input:: 
		prunedGraph: a Networkx Graph representing all the cloud clusters
		potentialMCCList: a list of dictionaries representing the possible MCCs within a path
		node: a string representing the cloud element currently being assessed
		CounterCriteriaAFlag: a boolean value indicating whether the node meets the MCC criteria A according to Laurent et al
		CounterCriteriaBFlag: a boolean value indicating whether the node meets the MCC criteria B according to Laurent et al
	
	Output:: 
		potentialMCCList: list of dictionaries representing all possible MCC within the path
			 dictionary = {"possMCCList":[(node,'I')], "fullMCSMCC":[(node,'I')], "CounterCriteriaA": CounterCriteriaA, "durationAandB": durationAandB}

	'''
	existingFrameFlag = False
	existingMCSFrameFlag = False
	predecessorsFlag = False
	predecessorsMCSFlag = False
	successorsFlag = False
	successorsMCSFlag = False
	frameNum = 0

	frameNum = int((node.split('CE')[0]).split('F')[1])
	if potentialMCCList==[]:
		#list empty
		stage = 'I'
		if CounterCriteriaAFlag == True and CounterCriteriaBFlag ==True:
			potentialMCCList.append({"possMCCList":[(node,stage)], "fullMCSMCC":[(node,stage)], "CounterCriteriaA": 1, "durationAandB": 1, "highestMCCnode":node, "frameNum":frameNum})	
		elif CounterCriteriaAFlag == True and CounterCriteriaBFlag == False:
			potentialMCCList.append({"possMCCList":[], "fullMCSMCC":[(node,stage)], "CounterCriteriaA": 1, "durationAandB": 0, "highestMCCnode":"", "frameNum":0})	
		elif CounterCriteriaAFlag == False and CounterCriteriaBFlag == False:
			potentialMCCList.append({"possMCCList":[], "fullMCSMCC":[(node,stage)], "CounterCriteriaA": 0, "durationAandB": 0, "highestMCCnode":"", "frameNum":0})	

	else:
		#list not empty
		predecessorsFlag, index = isThereALink(prunedGraph, 1,node,potentialMCCList,1)
		
		if predecessorsFlag == True:	

			for eachNode in potentialMCCList[index]["possMCCList"]:
				if int((eachNode[0].split('CE')[0]).split('F')[1]) == frameNum :
					existingFrameFlag = True
					
			#this MUST come after the check for the existing frame
			if CounterCriteriaAFlag == True and CounterCriteriaBFlag ==True:
				stage = 'M'
				potentialMCCList[index]["possMCCList"].append((node,stage))
				potentialMCCList[index]["fullMCSMCC"].append((node,stage))

			
			if existingFrameFlag == False:
				if CounterCriteriaAFlag == True and CounterCriteriaBFlag ==True:
					stage ='M'
					potentialMCCList[index]["CounterCriteriaA"]+= 1
					potentialMCCList[index]["durationAandB"]+=1
					if frameNum > potentialMCCList[index]["frameNum"]:
						potentialMCCList[index]["frameNum"] = frameNum
						potentialMCCList[index]["highestMCCnode"] = node
					return potentialMCCList

				#if this frameNum doesn't exist and this frameNum is less than the MCC node max frame Num (including 0), then append to fullMCSMCC list
				if frameNum > potentialMCCList[index]["frameNum"] or potentialMCCList[index]["frameNum"]==0:
					stage = 'I'
					if CounterCriteriaAFlag == True and CounterCriteriaBFlag == False:
						potentialMCCList.append({"possMCCList":[], "fullMCSMCC":[(node,stage)], "CounterCriteriaA": 1, "durationAandB": 0, "highestMCCnode":"", "frameNum":0})	
						return potentialMCCList
					elif CounterCriteriaAFlag == False and CounterCriteriaBFlag == False:
						potentialMCCList.append({"possMCCList":[], "fullMCSMCC":[(node,stage)], "CounterCriteriaA": 0, "durationAandB": 0, "highestMCCnode":"", "frameNum":0})	
						return potentialMCCList

			#if predecessor and this frame number already exist in the MCC list, add the current node to the fullMCSMCC list
			if existingFrameFlag == True:
				if CounterCriteriaAFlag == True and CounterCriteriaBFlag == False:
					potentialMCCList[index]["fullMCSMCC"].append((node,stage))
					potentialMCCList[index]["CounterCriteriaA"] +=1
					return potentialMCCList
				if CounterCriteriaAFlag == False:
					potentialMCCList[index]["fullMCSMCC"].append((node,stage))	
					return potentialMCCList	
				
		if predecessorsFlag == False:
			successorsFlag, index = isThereALink(prunedGraph, 2,node,potentialMCCList,2)
			
			if successorsFlag == True:
				for eachNode in potentialMCCList[index]["possMCCList"]: 
					if int((eachNode[0].split('CE')[0]).split('F')[1]) == frameNum:
						existingFrameFlag = True
						
				if CounterCriteriaAFlag == True and CounterCriteriaBFlag == True:
					stage = 'M'
					potentialMCCList[index]["possMCCList"].append((node,stage))
					potentialMCCList[index]["fullMCSMCC"].append((node,stage))
					if frameNum > potentialMCCList[index]["frameNum"] or potentialMCCList[index]["frameNum"] == 0:
						potentialMCCList[index]["frameNum"] = frameNum
						potentialMCCList[index]["highestMCCnode"] = node
					return potentialMCCList
		
				
				if existingFrameFlag == False:
					if stage == 'M':
						stage = 'D'
					if CounterCriteriaAFlag == True and CounterCriteriaBFlag ==True:
						potentialMCCList[index]["CounterCriteriaA"]+= 1
						potentialMCCList[index]["durationAandB"]+=1
					elif CounterCriteriaAFlag == True:
						potentialMCCList[index]["CounterCriteriaA"] += 1
					elif CounterCriteriaAFlag == False:
						potentialMCCList[index]["fullMCSMCC"].append((node,stage))
						return potentialMCCList
						#if predecessor and this frame number already exist in the MCC list, add the current node to the fullMCSMCC list
				else:
					if CounterCriteriaAFlag == True and CounterCriteriaBFlag == False:
						potentialMCCList[index]["fullMCSMCC"].append((node,stage))
						potentialMCCList[index]["CounterCriteriaA"] +=1
						return potentialMCCList
					if CounterCriteriaAFlag == False:
						potentialMCCList[index]["fullMCSMCC"].append((node,stage))	
						return potentialMCCList			

		#if this node isn't connected to exisiting MCCs check if it is connected to exisiting MCSs ...
		if predecessorsFlag == False and successorsFlag == False:
			stage = 'I'
			predecessorsMCSFlag, index = isThereALink(prunedGraph, 1,node,potentialMCCList,2)
			if predecessorsMCSFlag == True:
				if CounterCriteriaAFlag == True and CounterCriteriaBFlag == True:
					potentialMCCList[index]["possMCCList"].append((node,'M'))
					potentialMCCList[index]["fullMCSMCC"].append((node,'M'))
					potentialMCCList[index]["durationAandB"] += 1
					if frameNum > potentialMCCList[index]["frameNum"]:
						potentialMCCList[index]["frameNum"] = frameNum
						potentialMCCList[index]["highestMCCnode"] = node
					return potentialMCCList

				if potentialMCCList[index]["frameNum"] == 0 or frameNum <= potentialMCCList[index]["frameNum"]:
					if CounterCriteriaAFlag == True and CounterCriteriaBFlag == False:
						potentialMCCList[index]["fullMCSMCC"].append((node,stage))
						potentialMCCList[index]["CounterCriteriaA"] +=1
						return potentialMCCList
					elif CounterCriteriaAFlag == False:
						potentialMCCList[index]["fullMCSMCC"].append((node,stage))
						return potentialMCCList
			else:
				successorsMCSFlag, index = isThereALink(prunedGraph, 2,node,potentialMCCList,2)
				if successorsMCSFlag == True:
					if CounterCriteriaAFlag == True and CounterCriteriaBFlag == True:
						potentialMCCList[index]["possMCCList"].append((node,'M'))
						potentialMCCList[index]["fullMCSMCC"].append((node,'M'))
						potentialMCCList[index]["durationAandB"] += 1
						if frameNum > potentialMCCList[index]["frameNum"]:
							potentialMCCList[index]["frameNum"] = frameNum
							potentialMCCList[index]["highestMCCnode"] = node
						return potentialMCCList

					
					if potentialMCCList[index]["frameNum"] == 0 or frameNum <= potentialMCCList[index]["frameNum"]:
						if CounterCriteriaAFlag == True and CounterCriteriaBFlag == False:
							potentialMCCList[index]["fullMCSMCC"].append((node,stage))
							potentialMCCList[index]["CounterCriteriaA"] +=1
							return potentialMCCList
						elif CounterCriteriaAFlag == False:
							potentialMCCList[index]["fullMCSMCC"].append((node,stage))
							return potentialMCCList
					
			#if this node isn't connected to existing MCCs or MCSs, create a new one ...
			if predecessorsFlag == False and predecessorsMCSFlag == False and successorsFlag == False and successorsMCSFlag == False:	
				if CounterCriteriaAFlag == True and CounterCriteriaBFlag ==True:
					potentialMCCList.append({"possMCCList":[(node,stage)], "fullMCSMCC":[(node,stage)], "CounterCriteriaA": 1, "durationAandB": 1, "highestMCCnode":node, "frameNum":frameNum})	
				elif CounterCriteriaAFlag == True and CounterCriteriaBFlag == False:
					potentialMCCList.append({"possMCCList":[], "fullMCSMCC":[(node,stage)], "CounterCriteriaA": 1, "durationAandB": 0, "highestMCCnode":"", "frameNum":0})	
				elif CounterCriteriaAFlag == False and CounterCriteriaBFlag == False:
					potentialMCCList.append({"possMCCList":[], "fullMCSMCC":[(node,stage)], "CounterCriteriaA": 0, "durationAandB": 0, "highestMCCnode":"", "frameNum":0})	

	return potentialMCCList
#******************************************************************
def isThereALink(prunedGraph, upOrDown,node,potentialMCCList,whichList):
	'''
	Purpose:: 
		Utility script for updateMCCList mostly because there is no Pythonic way to break out of nested loops
	
	Input:: 
		prunedGraph:a Networkx Graph representing all the cloud clusters
		upOrDown: an integer representing 1- to do predecesor check and 2 - to do successor checkedNodesMCC
		node: a string representing the cloud element currently being assessed
		potentialMCCList: a list of dictionaries representing the possible MCCs within a path
		whichList: an integer representing which list ot check in the dictionary; 1- possMCCList, 2- fullMCSMCC
			
	Output:: 
		thisFlag: a boolean representing whether the list passed has in the parent or child of the node
		index: an integer representing the location in the potentialMCCList where thisFlag occurs

	'''
	thisFlag = False
	index = -1
	checkList =""
	if whichList == 1:
		checkList = "possMCCList"
	elif whichList ==2:
		checkList = "fullMCSMCC"

	#check parents
	if upOrDown == 1:
		for aNode in prunedGraph.predecessors(node):
			#reset the index counter for this node search through potentialMCCList
			index = -1
			for MCCDict in potentialMCCList:
				index += 1
				if aNode in list(x[0] for x in MCCDict[checkList]): 
					thisFlag = True
					#get out of looping so as to avoid the flag being written over when another node in the predecesor list is checked
					return thisFlag, index

	#check children
	if upOrDown == 2:
		for aNode in prunedGraph.successors(node):
			#reset the index counter for this node search through potentialMCCList
			index = -1
			for MCCDict in potentialMCCList:
				index += 1
				
				if aNode in list(x[0] for x in MCCDict[checkList]): 
					thisFlag = True
					return thisFlag, index

	return thisFlag, index
#******************************************************************
def maxExtentAndEccentricity(eachList):
	'''
	Purpose:: 
		Perform the final check for MCC based on maximum extent and eccentricity criteria

	Input:: 
		eachList: a list of strings  representing the node of the possible MCCs within a path

	Output:: 
		maxShieldNode: a string representing the node with the maximum maxShieldNode
	    definiteMCCFlag: a boolean indicating that the MCC has met all requirements

	'''
	maxShieldNode =''
	maxShieldArea = 0.0
	maxShieldEccentricity = 0.0
	definiteMCCFlag = False
	
	if eachList:
		for eachNode in eachList:
			if (thisDict(eachNode)['nodeMCSIdentifier'] == 'M' or thisDict(eachNode)['nodeMCSIdentifier'] == 'D') and thisDict(eachNode)['cloudElementArea'] > maxShieldArea:
				maxShieldNode = eachNode
				maxShieldArea = thisDict(eachNode)['cloudElementArea']
				
		maxShieldEccentricity = thisDict(maxShieldNode)['cloudElementEccentricity']
		if thisDict(maxShieldNode)['cloudElementEccentricity'] >= ECCENTRICITY_THRESHOLD_MIN and thisDict(maxShieldNode)['cloudElementEccentricity'] <= ECCENTRICITY_THRESHOLD_MAX :
			#criteria met
			definiteMCCFlag = True
			
	return maxShieldNode, definiteMCCFlag		
#******************************************************************
def findMaxDepthAndMinPath (thisPathDistanceAndLength):
	'''
	Purpose:: 
		To determine the maximum depth and min path for the headnode

	Input:: 
		tuple of dictionaries representing the shortest distance and paths for a node in the tree as returned by nx.single_source_dijkstra
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
			if pathLength == maxPathLength :
				if pathDistance <= minPath:
					minPath = pathLength
					#store details if absolute minPath and deepest
					minDistanceAndMaxPath = (pathDistance, path)
	return minDistanceAndMaxPath
#******************************************************************
def thisDict (thisNode):
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
#******************************************************************
def checkCriteria (thisCloudElementLatLon, aTemperature):
	'''
	Purpose:: 
		Determine if criteria B is met for a CEGraph

	Input:: 
		thisCloudElementLatLon: 2D array of (lat,lon) variable from the node dictionary being currently considered
		aTemperature:a integer representing the temperature maximum for masking

	Output :: 
		cloudElementArea: a floating-point number representing the area in the array that meet the criteria - criteriaB

	'''
	cloudElementCriteriaBLatLon=[]

	frame, CEcounter = ndimage.measurements.label(thisCloudElementLatLon, structure=STRUCTURING_ELEMENT)
	frameCEcounter=0
	#determine min and max values in lat and lon, then use this to generate teh array from LAT,LON meshgrid
	
	minLat = min(x[0] for x in thisCloudElementLatLon)
	maxLat = max(x[0]for x in thisCloudElementLatLon)
	minLon = min(x[1]for x in thisCloudElementLatLon)
	maxLon = max(x[1]for x in thisCloudElementLatLon)

	minLatIndex = np.argmax(LAT[:,0] == minLat)
	maxLatIndex = np.argmax(LAT[:,0]== maxLat)
	minLonIndex = np.argmax(LON[0,:] == minLon)
	maxLonIndex = np.argmax(LON[0,:] == maxLon)

	criteriaBframe = ma.zeros(((abs(maxLatIndex - minLatIndex)+1), (abs(maxLonIndex - minLonIndex)+1)))
	
	for x in thisCloudElementLatLon:
		#to store the values of the subset in the new array, remove the minLatIndex and minLonindex from the
		#index given in the original array to get the indices for the new array
		criteriaBframe[(np.argmax(LAT[:,0] == x[0]) - minLatIndex),(np.argmax(LON[0,:] == x[1]) - minLonIndex)] = x[2]

	#keep only those values < aTemperature
	tempMask = ma.masked_array(criteriaBframe, mask=(criteriaBframe >= aTemperature), fill_value = 0)
	
	#get the actual values that the mask returned
	criteriaB = ma.zeros((criteriaBframe.shape)).astype('int16')
	
	for index, value in utils.maenumerate(tempMask): 
		lat_index, lon_index = index			
		criteriaB[lat_index, lon_index]=value	

   	for count in xrange(CEcounter):
   		#[0] is time dimension. Determine the actual values from the data
   		#loc is a masked array
   		#***** returns elements down then across thus (6,4) is 6 arrays deep of size 4
   		try:

	   		loc = ndimage.find_objects(criteriaB)[0]
	   	except:
	   		#this would mean that no objects were found meeting criteria B
	   		print "no objects at this temperature!"
	   		cloudElementArea = 0.0
	   		return cloudElementArea, cloudElementCriteriaBLatLon
	   
	   	try:
	   		cloudElementCriteriaB = ma.zeros((criteriaB.shape))
	   		cloudElementCriteriaB =criteriaB[loc] 
	   	except:
	   		print "YIKESS"
	   		print "CEcounter ", CEcounter, criteriaB.shape
	   		print "criteriaB ", criteriaB

   		for index,value in np.ndenumerate(cloudElementCriteriaB):
   			if value !=0:
   				t,lat,lon = index
   				#add back on the minLatIndex and minLonIndex to find the true lat, lon values
   				lat_lon_tuple = (LAT[(lat),0], LON[0,(lon)],value)
   				cloudElementCriteriaBLatLon.append(lat_lon_tuple)

		cloudElementArea = np.count_nonzero(cloudElementCriteriaB)*XRES*YRES
		#do some cleaning up
		tempMask =[]
		criteriaB =[]
		cloudElementCriteriaB=[]

		return cloudElementArea, cloudElementCriteriaBLatLon
#******************************************************************
def hasMergesOrSplits (nodeList):
	'''
	Purpose:: 
		Determine if nodes within a path defined from shortest_path splittingNodeDict
	Input:: 
		nodeList: list of strings representing the nodes from a path
	Output:: 
		splitList: a list of strings representing all the nodes in the path that split
		mergeList: a list of strings representing all the nodes in the path that merged
	'''
	mergeList=[]
	splitList=[]

	for node,numParents in PRUNED_GRAPH.in_degree(nodeList).items():
		if numParents > 1:
			mergeList.append(node)

	for node, numChildren in PRUNED_GRAPH.out_degree(nodeList).items():
		if numChildren > 1:
			splitList.append(node)
	#sort
	splitList.sort(key=lambda item:(len(item.split('C')[0]), item.split('C')[0]))
	mergeList.sort(key=lambda item:(len(item.split('C')[0]), item.split('C')[0]))
			
	return mergeList,splitList
#******************************************************************
def allAncestors(path, aNode):
	'''
	Purpose:: 
		Utility script to provide the path leading up to a nodeList

	Input:: 
		path: a list of strings representing the nodes in the path 
	    aNode: a string representing a node to be checked for parents

	Output:: 
		path: a list of strings representing the list of the nodes connected to aNode through its parents
		numOfChildren: an integer representing the number of parents of the node passed
	'''

	numOfParents = PRUNED_GRAPH.in_degree(aNode)
	try:
		if PRUNED_GRAPH.predecessors(aNode) and numOfParents <= 1:
			path = path + PRUNED_GRAPH.predecessors(aNode)
			thisNode = PRUNED_GRAPH.predecessors(aNode)[0]
			return allAncestors(path,thisNode)
		else:
			path = path+aNode
			return path, numOfParents
	except:
		return path, numOfParents
#******************************************************************
def allDescendants(path, aNode):
	'''
	Purpose:: 
		Utility script to provide the path leading up to a nodeList

	Input:: 
		path: a list of strings representing the nodes in the path 
	    aNode: a string representing a node to be checked for children

	Output:: 
		path: a list of strings representing the list of the nodes connected to aNode through its children
		numOfChildren: an integer representing the number of children of the node passed
	'''

	numOfChildren = PRUNED_GRAPH.out_degree(aNode)
	try:
		if PRUNED_GRAPH.successors(aNode) and numOfChildren <= 1:
			path = path + PRUNED_GRAPH.successors(aNode)
			thisNode = PRUNED_GRAPH.successors(aNode)[0]
			return allDescendants(path,thisNode)
		else:
			path = path + aNode
			#i.e. PRUNED_GRAPH.predecessors(aNode) is empty
			return path, numOfChildren
	except:
		#i.e. PRUNED_GRAPH.predecessors(aNode) threw an exception
		return path, numOfChildren
#******************************************************************
def addInfothisDict (thisNode, cloudElementArea,criteriaB):
	'''
	Purpose:: 
		Update original dictionary node with information

	Input:: 
		thisNode: a string representing the unique ID of a node
		cloudElementArea: a floating-point number representing the area of the cloud element
		criteriaB: a masked array of floating-point numbers representing the lat,lons meeting the criteria  

	Output:: None 

	'''
	for eachdict in CLOUD_ELEMENT_GRAPH.nodes(thisNode):
		if eachdict[1]['uniqueID'] == thisNode:
			eachdict[1]['CriteriaBArea'] = cloudElementArea
			eachdict[1]['CriteriaBLatLon'] = criteriaB
	return
#******************************************************************
def addNodeBehaviorIdentifier (thisNode, nodeBehaviorIdentifier):
	'''
	Purpose:: add an identifier to the node dictionary to indicate splitting, merging or neither node

	Input:: 
		thisNode: a string representing the unique ID of a node
	    nodeBehaviorIdentifier: a string representing the behavior S- split, M- merge, B- both split and merge, N- neither split or merge 

	Output :: None

	'''
	for eachdict in CLOUD_ELEMENT_GRAPH.nodes(thisNode):
		if eachdict[1]['uniqueID'] == thisNode:
			if not 'nodeBehaviorIdentifier' in eachdict[1].keys():
				eachdict[1]['nodeBehaviorIdentifier'] = nodeBehaviorIdentifier
	return
#******************************************************************
def addNodeMCSIdentifier (thisNode, nodeMCSIdentifier):
	'''
	Purpose:: 
		Add an identifier to the node dictionary to indicate splitting, merging or neither node

	Input:: 
		thisNode: a string representing the unique ID of a node
		nodeMCSIdentifier: a string representing the stage of the MCS lifecyle  'I' for Initiation, 'M' for Maturity, 'D' for Decay

	Output :: None

	'''
	for eachdict in CLOUD_ELEMENT_GRAPH.nodes(thisNode):
		if eachdict[1]['uniqueID'] == thisNode:
			if not 'nodeMCSIdentifier' in eachdict[1].keys():
				eachdict[1]['nodeMCSIdentifier'] = nodeMCSIdentifier
	return
#******************************************************************
def updateNodeMCSIdentifier (thisNode, nodeMCSIdentifier):
	'''
	Purpose:: 
		Update an identifier to the node dictionary to indicate splitting, merging or neither node

	Input:: 
		thisNode: thisNode: a string representing the unique ID of a node
		nodeMCSIdentifier: a string representing the stage of the MCS lifecyle  'I' for Initiation, 'M' for Maturity, 'D' for Decay  

	Output :: None

	'''
	for eachdict in CLOUD_ELEMENT_GRAPH.nodes(thisNode):
		if eachdict[1]['uniqueID'] == thisNode:
			eachdict[1]['nodeMCSIdentifier'] = nodeBehaviorIdentifier

	return
#******************************************************************
def eccentricity (cloudElementLatLon):
	'''
	Purpose::
	    Determines the eccentricity (shape) of contiguous boxes 
	    Values tending to 1 are more circular by definition, whereas 
	    values tending to 0 are more linear
	
	Input::
		cloudElementLatLon: 2D array in (lat,lon) representing T_bb contiguous squares 
		
	Output::
		epsilon: a floating-point representing the eccentricity of the matrix passed
	
	'''
	
	epsilon = 0.0
	
	#loop over all lons and determine longest (non-zero) col
	#loop over all lats and determine longest (non-zero) row
	for latLon in cloudElementLatLon:
	    #assign a matrix to determine the legit values
	    
	    nonEmptyLons = sum(sum(cloudElementLatLon)>0)
        nonEmptyLats = sum(sum(cloudElementLatLon.transpose())>0)
        
        lonEigenvalues = 1.0 * nonEmptyLats / (nonEmptyLons+0.001) #for long oval on y axis
        latEigenvalues = 1.0 * nonEmptyLons / (nonEmptyLats +0.001) #for long oval on x-axs
        epsilon = min(latEigenvalues,lonEigenvalues)
        
	return epsilon
#******************************************************************
def cloudElementOverlap (currentCELatLons, previousCELatLons):
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

	latlonprev =[]
	latloncurr = []
	count = 0 
	percentageOverlap = 0.0
	areaOverlap = 0.0

	#remove the temperature from the tuples for currentCELatLons and previousCELatLons then check for overlap
	latlonprev = [(x[0],x[1]) for x in previousCELatLons]
	latloncurr = [(x[0],x[1]) for x in currentCELatLons]  

	#find overlap
	count = len(list(set(latloncurr)&set(latlonprev)))

	#find area overlap
	areaOverlap = count*XRES*YRES
	
	#find percentage
	percentageOverlap = max(((count*1.0)/(len(latloncurr)*1.0)),((count*1.0)/(len(latlonprev)*1.0)))
	
	return percentageOverlap, areaOverlap
#******************************************************************
