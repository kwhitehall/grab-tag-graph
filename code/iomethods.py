import datetime
from datetime import timedelta, datetime
import glob
import itertools
from netCDF4 import Dataset, num2date, date2num
import numpy as np
import numpy.ma as ma
import os

import utils
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
def readMergData(dirname, filelist=None):
	'''
	Purpose::
	    Read MERG data into RCMES format
	
	Input::
	    dirname: a string representing the directory to the MERG files in NETCDF format
	    filelist (optional): a list of strings representing the filenames betweent the start and end dates provided
	
	Output::
	    A 3D masked array (t,lat,lon) with only the variables which meet the minimum temperature 
	    criteria for each frame

	Assumptions::
	    The MERG data has been converted to NETCDF using LATS4D
	    The data has the same lat/lon format

	TODO:: figure out how to use netCDF4 to do the clipping tmp = netCDF4.Dataset(filelist[0])

	'''

	global LAT
	global LON

	# these strings are specific to the MERG data
	mergVarName = 'ch4'
	mergTimeVarName = 'time'
	mergLatVarName = 'latitude'
	mergLonVarName = 'longitude'
	
	filelistInstructions = dirname + '/*'
	#filelist = glob.glob(filelistInstructions)
	if filelist == None:
		filelist = glob.glob(filelistInstructions)

	
	#sat_img is the array that will contain all the masked frames
	mergImgs = []
	#timelist of python time strings
	timelist = [] 
	time2store = None
	tempMaskedValueNp =[]
	

	filelist.sort()
	nfiles = len(filelist)

	# Crash nicely if there are no netcdf files
	if nfiles == 0:
		print 'Error: no files in this directory! Exiting elegantly'
		sys.exit()
	else:
		# Open the first file in the list to read in lats, lons and generate the  grid for comparison
		#tmp = Nio.open_file(filelist[0], format='nc')
		tmp = Dataset(filelist[0], 'r+',format='NETCDF4')

		alllatsraw = tmp.variables[mergLatVarName][:]
		alllonsraw = tmp.variables[mergLonVarName][:]
		alllonsraw[alllonsraw > 180] = alllonsraw[alllonsraw > 180] - 360.  # convert to -180,180 if necessary
		
		#get the lat/lon info data (different resolution)
		latminNETCDF = utils.findNearest(alllatsraw, float(LATMIN))
		latmaxNETCDF = utils.findNearest(alllatsraw, float(LATMAX))
		lonminNETCDF = utils.findNearest(alllonsraw, float(LONMIN))
		lonmaxNETCDF = utils.findNearest(alllonsraw, float(LONMAX))
		latminIndex = (np.where(alllatsraw == latminNETCDF))[0][0]
		latmaxIndex = (np.where(alllatsraw == latmaxNETCDF))[0][0]
		lonminIndex = (np.where(alllonsraw == lonminNETCDF))[0][0]
		lonmaxIndex = (np.where(alllonsraw == lonmaxNETCDF))[0][0]
		
		#subsetting the data
		latsraw = alllatsraw[latminIndex: latmaxIndex]
		lonsraw = alllonsraw[lonminIndex:lonmaxIndex]

		LON, LAT = np.meshgrid(lonsraw, latsraw)
		#clean up
		latsraw =[]
		lonsraw = []
		nygrd = len(LAT[:, 0]); nxgrd = len(LON[0, :])
		tmp.close
	
	for files in filelist:
		
		try:
			thisFile = Dataset(files,'r', format='NETCDF4')
			#clip the dataset according to user lat,lon coordinates
			#mask the data and fill with zeros for later 
			tempRaw = thisFile.variables[mergVarName][:,latminIndex:latmaxIndex,lonminIndex:lonmaxIndex].astype('int16')
 			tempMask = ma.masked_array(tempRaw, mask=(tempRaw > T_BB_MAX), fill_value=0) 
 			#get the actual values that the mask returned
			tempMaskedValue = ma.zeros((tempRaw.shape)).astype('int16')
			
			for index, value in utils.maenumerate(tempMask): 
				time_index, lat_index, lon_index = index		
				tempMaskedValue[time_index,lat_index, lon_index]=value	
			
			xtimes = thisFile.variables[mergTimeVarName]
			#convert this time to a python datastring
			time2store, _ = utils.getModelTimes(xtimes, mergTimeVarName)
			#extend instead of append because getModelTimes returns a list already and we don't 
			#want a list of list
			timelist.extend(time2store)
			mergImgs.extend(tempMaskedValue) 
			thisFile.close
			thisFile = None
			
		except:
			print "bad file! ", files

	mergImgs = ma.array(mergImgs)
	
	return mergImgs, timelist, LAT, LON
#******************************************************************
