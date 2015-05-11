
import glob
import numpy.ma as ma
import numpy as np
import re
import string
import os

from netCDF4 import Dataset
from datetime import timedelta, datetime
from os import path

import utils
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
STRUCTURING_ELEMENT = [[0,1,0],[1,1,1],[0,1,0]] #the matrix for determining the pattern for the contiguous boxes and must
                                                #have same rank of the matrix it is being compared against
#criteria for determining cloud elements and edges
T_BB_MAX = 243  #warmest temp to allow (-30C to -55C according to Morel and Sensi 2002)
T_BB_MIN = 218  #cooler temp for the center of the system
CONVECTIVE_FRACTION = 0.90 #the min temp/max temp that would be expected in a CE.. this is highly conservative (only a 10K difference)
MIN_MCS_DURATION = 3    #minimum time for a MCS to exist
AREA_MIN = 2400.0       #minimum area for CE criteria in km^2 according to Vila et al. (2008) is 2400
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


def checkForFiles(dirPath, startTime, endTime, tdelta, tRes):
    '''
        Purpose:: To ensure all the files between the starttime and endTime
                  exist in the directory supplied

        Input:: dirPath: a string representing the path to the files 

                startTime: a string representing the starttime.
                    Must have atleast yyyymm
                
                endTime: a string representing the endTime
                    Must have atleast yyyymm

                tdelta: an integer representing the time between files

                tRes: a string representing the time resolution for tdelta. 
                    Must be month, day, hour or minutes

        Output::
                status: a boolean representing whether all files exists

        Assumptions:: the filename contains the time

    '''
    filelist = []
    filenamePattern = ''
    startFile = ''
    endFile = ''
    currFile = ''
    hasDelimiter = False
    startTimeInFile = ''
    endTimeInFile = ''
    currTimeInfile = ''
    tokenCounter = 0
    
    if 'month' in tRes:
        currFileTime = datetime.strptime(startTime[:6], "%Y%m")
        tRes = 'month'
    if 'day' in tRes:
        currFileTime = datetime.strptime(startTime[:8], "%Y%m%d")
        tRes = "day"
    if 'hour' in tRes:
        currFileTime = datetime.strptime(startTime[:10], "%Y%m%d%H")
        tRes = 'hour'
    if 'minute' in tRes:
        currFileTime = datetime.strptime(startTime[:12], "%Y%m%d%H%m")
        tRes = "minute"

    
    filelist = filter(path.isfile, glob.glob((dirPath+'/*.nc')))
    filelist.sort()
    
    #check for the filename pattern 
    for eachPart in re.split(r'[_,-,.,/]',re.split(r'.nc',path.basename(filelist[0]))[0]):
        tokenCounter += 1
        if tokenCounter == 1:
            filenamePattern += eachPart
        if eachPart.isdigit():
            if len(eachPart) >= 6:
                hasDelimiter = True
                startTimeInFile += eachPart +'*'
            elif eachPart in startTime:
                startTimeInFile += eachPart +'*'

    #can glob also match a file pattern? i.e. in this case, the date needs to have a delimiter in order to be matched

    if hasDelimiter is False:
        fileDate = int(re.search(r'\d+',re.split(r'.nc',path.basename(filelist[0]))[0]).group())
        filenamePattern = re.split(str(fileDate),path.basename(filelist[0]))[0]

    startFile = glob.glob(dirPath+'/'+filenamePattern +'*'+startTimeInFile)[0]
    endTimeInFile = find_time_in_file(endTime, startTimeInFile)
    endFile = glob.glob(dirPath+'/'+filenamePattern + '*'+endTimeInFile+'*')[0]
    
    print "filenamePattern is ", filenamePattern, startTime, startTimeInFile
    print "dirPath is ", dirPath
    print "startFile is ", startFile
    print "endFile is ", endFile

    currFile = startFile
    filelist =[]
    
    #check for files between startTime and endTime
    while currFile is not endFile:
        if not path.isfile(currFile):
            status = False
            return status, filelist
        else:
            filelist.append(currFile)
    
        status = True
        if currFile == endFile:
            break

        if 'month' in tRes:
            currFileTime += timedelta(days = 31*tdelta)
            currTimeInFile = find_time_in_file(currFileTime.strftime("%Y%m"), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]
        if 'day' in tRes:
            currFileTime += timedelta(days = tdelta)
            currTimeInFile = find_time_in_file(currFileTime.strftime("%Y%m%d"), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]
        if 'hour' in tRes:
            currFileTime += timedelta(hours=tdelta)
            currTimeInFile = find_time_in_file(currFileTime.strftime("%Y%m%d%H"), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]
        if 'minute' in tRes:
            currFileTime += timedelta(minutes=tdelta)
            currTimeInFile = find_time_in_file(currFileTime.strftime("%Y%m%d%H%M"), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]

    return status,filelist
#******************************************************************
def find_time_in_file(myTime, myTimeInFile):
    '''
        Purpose:: To return the file pattern of the time string 

        Inputs:: myTime: a string in time format representing the time 
                 myTimeInFile: a string representing the time pattern of the time 
                 in a file

        Returns:: currTimeInfile: a string representing the pattern of the time in the file

    '''

    lastPos = 0
    currTimeInFile = ''
    
    for eachPart in re.split(r'[*]',myTimeInFile):
        if eachPart:
            currTimeInFile += myTime[ lastPos:lastPos+len(eachPart) ]+'*'
            lastPos += len(eachPart)
            
    return currTimeInFile
#******************************************************************
def createMainDirectory(mainDirStr):
    '''
        Purpose:: To create the main directory for storing information and
            the subdirectories for storing information
        
        Input:: mainDir: string representing the directory for where all 
                information generated from the program are to be stored
        
        Returns:: None

        Outputs:: a file structure where data generated from GTG will be stored 

        Assumptions:: The user running the program can write at mainDirStr 

    '''
    global MAINDIRECTORY

    MAINDIRECTORY = mainDirStr
    #if directory doesnt exist, create it
    if not os.path.exists(MAINDIRECTORY):
        os.makedirs(MAINDIRECTORY)

    os.chdir((MAINDIRECTORY))
    #create the subdirectories
    try:
        os.makedirs('images')
        os.makedirs('textFiles')
        os.makedirs('MERGnetcdfCEs')
        os.makedirs('TRMMnetcdfCEs')
    except:
        print "Directory exists already!!!"
        #TODO: some nice way of prompting if it is ok to continue...or just leave

    return MAINDIRECTORY
#******************************************************************
def readData(dirname, varName, latName, lonName, filelist=None):
    '''
        Purpose::
            Read gridded data into (t, lat, lon) arrays for processing

        Inputs::
            dirname: a string representing the directory to the MERG files in NETCDF format
            varName: a string representing the variable name to use from the file
            latName: a string representing the latitude from the file's metadata
            lonName: a string representing the longitude from the file's metadata
            filelist (optional): a list of strings representing the filenames betweent the start and end dates provided

        Returns:

        Outputs::
            A 3D masked array (t,lat,lon) with only the variables which meet the minimum temperature
            criteria for each frame

        Assumptions::
            (1) All the files requested to extract data are from the same instrument/model, and thus have the same metadata
            properties (varName, latName, lonName) as entered 
            (2) Assumes rectilinear grids for input datasets i.e. lat, lon will be 1D arrays
    '''

    global LAT
    global LON

    timeName = 'time'

    filelistInstructions = dirname + '/*'
    if filelist == None:
        filelist = glob.glob(filelistInstructions)


    inputData = []
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
        tmp = Dataset(filelist[0], 'r+',format='NETCDF4')

        alllatsraw = tmp.variables[latName][:]
        alllonsraw = tmp.variables[lonName][:]
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
        
        latsraw =[]
        lonsraw = []
        nygrd = len(LAT[:, 0]); nxgrd = len(LON[0, :])
        tmp.close

    for files in filelist:

        try:
            thisFile = Dataset(files,'r', format='NETCDF4')
            #clip the dataset according to user lat,lon coordinates
            #mask the data and fill with zeros for later
            tempRaw = thisFile.variables[varName][:,latminIndex:latmaxIndex,lonminIndex:lonmaxIndex].astype('int16')
            tempMask = ma.masked_array(tempRaw, mask=(tempRaw > T_BB_MAX), fill_value=0)
            #get the actual values that the mask returned
            tempMaskedValue = ma.zeros((tempRaw.shape)).astype('int16')

            for index, value in utils.maenumerate(tempMask):
                time_index, lat_index, lon_index = index
                tempMaskedValue[time_index,lat_index, lon_index]=value

            xtimes = thisFile.variables[timeName]

            #convert this time to a python datastring
            time2store, _ = getModelTimes(xtimes, timeName)

            #extend instead of append because getModelTimes returns a list already and we don't
            #want a list of list
            timelist.extend(time2store)
            inputData.extend(tempMaskedValue)
            thisFile.close
            thisFile = None

        except:
            print "bad file! ", files

    inputData = ma.array(inputData)

    return inputData, timelist, LAT, LON
#******************************************************************
def getModelTimes(xtimes, timeVarName):
    '''
    Purpose:: Routine to convert from model times ('hours since 1900...', 'days since ...')
    into a python datetime structure. Leveraged from Apache OCW

    Inputs::
        modelFile: a string representing the path to the model tile you want to 
        extract the times list and modelTimeStep from
        timeVarName: a string representing the name of the time variable in the model file

    Returns::
        times: a list of python datetime objects describing model data times
        modelTimeStep: a string representing the time step found in the file e.g. 
        'hourly','daily','monthly','annual'

    Outputs:: None

    Assumptions:: None
    
    '''

    timeFormat = xtimes.units

    try:
        sinceLoc = re.search('since', timeFormat).end()

    except AttributeError:
        print 'Error decoding model times: time variable attributes do not contain "since"'
        raise

    units = None
    TIME_UNITS = ('minutes', 'hours', 'days', 'months', 'years')
    for unit in TIME_UNITS:
        if re.search(unit, timeFormat):
            units = unit
            break

    base_time_string = string.lstrip(timeFormat[sinceLoc:])
    base_time = decodeTimeFromString(base_time_string)

    times = []

    for xtime in xtimes[:]:

        #TODO: This casting may cause problems for data that is hourly with more than one timestep in it
        xtime = int(xtime)

        if xtime%10 != 0:
            xtime = 1

        if units == 'minutes':
            dt = timedelta(minutes=xtime)
            new_time = base_time + dt
        elif units == 'hours':
            dt = timedelta(hours=int(xtime))
            new_time = base_time + dt# timedelta(hours=int(xtime))
        elif units == 'days':
            dt = timedelta(days=xtime)
            new_time = base_time + dt
        elif units == 'months':
            # NB. adding months in python is complicated as month length varies and hence ambiguous.+
            # Perform date arithmetic manually
            #  Assumption: the base_date will usually be the first of the month
            #              NB. this method will fail if the base time is on the 29th or higher day of month
            #                      -as can't have, e.g. Feb 31st.
            new_month = int(base_time.month + xtime % 12)
            new_year = int(math.floor(base_time.year + xtime / 12.))
            new_time = datetime.datetime(new_year, new_month, base_time.day, base_time.hour, base_time.second, 0)
        elif units == 'years':
            dt = datetime.timedelta(years=xtime)
            new_time = base_time + dt

        times.append(new_time)

    try:
        if len(xtimes) == 1:
            timeStepLength = 0
        else:
            timeStepLength = int(xtimes[1] - xtimes[0] + 1.e-12)

        modelTimeStep = getModelTimeStep(units, timeStepLength)

    except:
        raise

    return times, modelTimeStep
#******************************************************************
def getModelTimeStep(units, stepSize):
    '''
        Purpose:: To determine the time intervals of input data. 
                  Levaraged from Apache OCW

        Inputs:: units: a string representing the time units found in the
                 file metadata
                 stepSize: an integer representing the time interval found
                 in the file's metadata

        Returns:: modelTimeStep: a string representing the step interval in the 
                 dataset e.g. 'hourly', 'daily', etc.
    '''

    if units == 'minutes':
        if stepSize == 60:
            modelTimeStep = 'hourly'
        elif stepSize == 1440:
            modelTimeStep = 'daily'
        # 28 days through 31 days
        elif 40320 <= stepSize <= 44640:
            modelTimeStep = 'monthly'
        # 365 days through 366 days
        elif 525600 <= stepSize <= 527040:
            modelTimeStep = 'annual'
        else:
            raise Exception('model data time step interval exceeds the max time interval (annual)', units, stepSize)

    elif units == 'hours':
        #TODO: need a check for fractional hrs and only one hr i.e. stepSize=0 e.g. with MERG data
        if stepSize == 0 or stepSize == 1:
            modelTimeStep = 'hourly'
        elif stepSize == 24:
            modelTimeStep = 'daily'
        elif 672 <= stepSize <= 744:
            modelTimeStep = 'monthly'
        elif 8760 <= stepSize <= 8784:
            modelTimeStep = 'annual'
        else:
            raise Exception('model data time step interval exceeds the max time interval (annual)', units, stepSize)

    elif units == 'days':
        if stepSize == 1:
            modelTimeStep = 'daily'
        elif 28 <= stepSize <= 31:
            modelTimeStep = 'monthly'
        elif 365 <= stepSize <= 366:
            modelTimeStep = 'annual'
        else:
            raise Exception('model data time step interval exceeds the max time interval (annual)', units, stepSize)

    elif units == 'months':
        if stepSize == 1:
            modelTimeStep = 'monthly'
        elif stepSize == 12:
            modelTimeStep = 'annual'
        else:
            raise Exception('model data time step interval exceeds the max time interval (annual)', units, stepSize)

    elif units == 'years':
        if stepSize == 1:
            modelTimeStep = 'annual'
        else:
            raise Exception('model data time step interval exceeds the max time interval (annual)', units, stepSize)

    else:
        errorMessage = 'the time unit ', units, ' is not currently handled in this version.'
        raise Exception(errorMessage)

    return modelTimeStep
#******************************************************************
def decodeTimeFromString(timeString):
    '''
       Purpose:: Decodes string into a python datetime object

       Inputs:: timeString: a string representing a date/time

       Returns:: myTime: a python datetime object of the time_string

    '''
    
    # This will deal with times that use decimal seconds
   if '.' in timeString:
        timeString = timeString.split('.')[0] + '0'

    else:
        pass

    try:
        if datetime.strptime(timeString,'%Y-%m-%d %H'):
            myTime = datetime.strptime(timeString,'%Y-%m-%d %H')
        elif datetime.strptime(timeString,'%Y-%m-%d %H%M'):
            myTime = datetime.strptime(timeString,'%Y-%m-%d %H%M')
        return myTime

    except ValueError:
        pass

    print 'Error decoding time string: string does not match a predefined time format'
    return 0
#******************************************************************

