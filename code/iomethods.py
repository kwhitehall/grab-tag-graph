import glob
import os
import re
import string
import subprocess
import sys
from datetime import timedelta, datetime
from os import path

import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

import utils
import variables


def check_for_files(dirPath, startTime, endTime, tdelta, tRes):
    '''
        Purpose:: To ensure all the files between the startTime and endTime
                  exist in the directory supplied

        Input:: dirPath: a string representing the path to the files

                startTime: a string representing the startTime.
                    Must have at least yyyymm

                endTime: a string representing the endTime
                    Must have at least yyyymm

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
    currTimeInFile = ''
    tokenCounter = 0

    if 'month' in tRes:
        currFileTime = datetime.strptime(startTime[:6], '%Y%m')
        tRes = 'month'
    if 'day' in tRes:
        currFileTime = datetime.strptime(startTime[:8], '%Y%m%d')
        tRes = 'day'
    if 'hour' in tRes:
        currFileTime = datetime.strptime(startTime[:10], '%Y%m%d%H')
        tRes = 'hour'
    if 'minute' in tRes:
        currFileTime = datetime.strptime(startTime[:12], '%Y%m%d%H%m')
        tRes = 'minute'

    filelist = filter(path.isfile, glob.glob((dirPath+'/*.nc')))
    filelist.sort()

    # check for the filename pattern
    for eachPart in re.split(r'[_,-,.,/]', re.split(r'.nc', path.basename(filelist[0]))[0]):
        tokenCounter += 1
        if tokenCounter == 1:
            filenamePattern += eachPart
        if eachPart.isdigit():
            if len(eachPart) >= 6:
                hasDelimiter = True
                startTimeInFile += eachPart + '*'
            elif eachPart in startTime:
                startTimeInFile += eachPart + '*'

    if hasDelimiter is False:
        fileDate = int(re.search(r'\d+', re.split(r'.nc', path.basename(filelist[0]))[0]).group())
        filenamePattern = re.split(str(fileDate), path.basename(filelist[0]))[0]

    startFile = glob.glob(dirPath+'/'+filenamePattern + '*' + startTimeInFile)[0]
    endTimeInFile = find_time_in_file(endTime, startTimeInFile)
    # print("filenamePattern: "+filenamePattern)
    # print("endTimeInFile: "+endTimeInFile)
    # print(dirPath+'/'+filenamePattern + '*'+endTimeInFile+'*')
    endFile = glob.glob(dirPath+'/'+filenamePattern + '*'+endTimeInFile+'*')[0]

    currFile = startFile
    filelist = []

    # Check for files between startTime and endTime
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
            currFileTime += timedelta(days=31*tdelta)
            currTimeInFile = find_time_in_file(currFileTime.strftime('%Y%m'), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]
        if 'day' in tRes:
            currFileTime += timedelta(days=tdelta)
            currTimeInFile = find_time_in_file(currFileTime.strftime('%Y%m%d'), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]
        if 'hour' in tRes:
            currFileTime += timedelta(hours=tdelta)
            currTimeInFile = find_time_in_file(currFileTime.strftime('%Y%m%d%H'), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]
        if 'minute' in tRes:
            currFileTime += timedelta(minutes=tdelta)
            currTimeInFile = find_time_in_file(currFileTime.strftime('%Y%m%d%H%M'), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]

    return status, filelist
# **********************************************************************************************************************
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

    for eachPart in re.split(r'[*]', myTimeInFile):
        if eachPart:
            currTimeInFile += myTime[lastPos:lastPos+len(eachPart)]+'*'
            lastPos += len(eachPart)

    return currTimeInFile
# ***********************************************************************************************************************
def read_vars(userVariables):
    # for GrADs
    subprocess.call('export DISPLAY=:0.0', shell=True)

    graphVariables = variables.define_graph_variables()
    # ---------------------------------- end user inputs --------------------------------------
    # Checks that inputs are ok

    try:
        if not os.path.exists(userVariables.DIRS['TRMMdirName']):
            print "Error: TRMM invalid path!"
            userVariables.DIRS['TRMMdirName'] = raw_input("> Please enter the location to the raw TRMM netCDF files: \n")
    except:
        pass

    try:
        if not os.path.exists(userVariables.DIRS['CEoriDirName']):
            print "Error! MERG invalid path!"
            userVariables.DIRS['CEoriDirName'] = raw_input("> Please enter the directory to the MERG netCDF files: \n")
    except:
        print "..."   

    # check validity of time
    while utils.valid_date(userVariables.startDateTime) != True:
        print "Invalid time entered for startDateTime!"

    while utils.valid_date(userVariables.endDateTime) != True:
        print "Invalid time entered for endDateTime!"
        
    # check if all the files exisits in the MERG and TRMM directories entered
    test, _ = check_for_files(userVariables.DIRS['TRMMdirName'], userVariables.startDateTime, userVariables.endDateTime, 3, 'hour')
    if test is False:
        print "Error with files in the TRMM directory entered. Please check your files before restarting. "
        return
    
    test, userVariables.filelist = check_for_files(userVariables.DIRS['CEoriDirName'], userVariables.startDateTime, userVariables.endDateTime, 1, 'hour')

    if test is False:
        print "Error with files in the original MERG directory entered. Please check your files before restarting. "
        return

    # create main directory and file structure for storing intel
    userVariables.DIRS['mainDirStr'] = create_main_directory(userVariables.DIRS['mainDirStr'])
    TRMMCEdirName = userVariables.DIRS['mainDirStr']+'/TRMMnetcdfCEs'
    CEdirName = userVariables.DIRS['mainDirStr']+'/MERGnetcdfCEs'

    return graphVariables
# **********************************************************************************************************************
def create_main_directory(mainDirStr):
    '''
        Purpose:: To create the main directory for storing information and
            the subdirectories for storing information

        Input:: mainDir: string representing the directory for where all
                information generated from the program are to be stored

        Returns:: None

        Outputs:: a file structure where data generated from GTG will be stored

        Assumptions:: The user running the program can write at mainDirStr

    '''
    global MAIN_DIRECTORY

    MAIN_DIRECTORY = mainDirStr
    # if directory doesnt exist, create it
    if not os.path.exists(MAIN_DIRECTORY):
        os.makedirs(MAIN_DIRECTORY)

    os.chdir((MAIN_DIRECTORY))
    # create the subdirectories
    try:
        os.makedirs('images')
        os.makedirs('textFiles')
        os.makedirs('MERGnetcdfCEs')
        os.makedirs('TRMMnetcdfCEs')
    except:
        print 'Directory exists already!!!'
        # TODO: some nice way of prompting if it is ok to continue...or just leave

    return MAIN_DIRECTORY
# **********************************************************************************************************************
def read_data(varName, latName, lonName, userVariables, filelist=None):
    '''
        Purpose::
            Read gridded data into (t, lat, lon) arrays for processing

        Inputs::
            varName: a string representing the variable name to use from the file
            latName: a string representing the latitude from the file's metadata
            lonName: a string representing the longitude from the file's metadata
            userVariables: 
            filelist (optional): a list of strings representing the filenames between the start and end dates provided

        Returns:

        Outputs::
            A 3D masked array (t,lat,lon) with only the variables which meet the minimum temperature
            criteria for each frame

        Assumptions::
            (1) All the files requested to extract data are from the same instrument/model, and thus have the same
            metadata properties (varName, latName, lonName) as entered
            (2) Assumes rectilinear grids for input datasets i.e. lat, lon will be 1D arrays
    '''

    global LAT
    global LON

    timeName = 'time'

    filelistInstructions = userVariables.DIRS['CEoriDirName']+'/*'
    
    if filelist is None and userVariables.filelist is None:
        userVariables.filelist = glob.glob(filelistInstructions)
        
    userVariables.filelist.sort()

    inputData = []
    timelist = []
    time2store = None
    tempMaskedValueNp = []

    nfiles = len(userVariables.filelist)

    # Crash nicely if there are no netcdf files
    if nfiles == 0:
        print 'Error: no files in this directory! Exiting elegantly'
        sys.exit()
    else:
        # Open the first file in the list to read in lats, lons and generate the  grid for comparison
        tmp = Dataset(userVariables.filelist[0], 'r+', format='NETCDF4')

        alllatsraw = tmp.variables[latName][:]
        alllonsraw = tmp.variables[lonName][:]
        alllonsraw[alllonsraw > 180] = alllonsraw[alllonsraw > 180] - 360.  # convert to -180,180 if necessary

        # Get the lat/lon info data (different resolution)
        latminNETCDF = utils.find_nearest(alllatsraw, float(userVariables.LATMIN))
        latmaxNETCDF = utils.find_nearest(alllatsraw, float(userVariables.LATMAX))
        lonminNETCDF = utils.find_nearest(alllonsraw, float(userVariables.LONMIN))
        lonmaxNETCDF = utils.find_nearest(alllonsraw, float(userVariables.LONMAX))
        latminIndex = (np.where(alllatsraw == latminNETCDF))[0][0]
        latmaxIndex = (np.where(alllatsraw == latmaxNETCDF))[0][0]
        lonminIndex = (np.where(alllonsraw == lonminNETCDF))[0][0]
        lonmaxIndex = (np.where(alllonsraw == lonmaxNETCDF))[0][0]

        # Subsetting the data
        latsraw = alllatsraw[latminIndex: latmaxIndex]
        lonsraw = alllonsraw[lonminIndex:lonmaxIndex]

        LON, LAT = np.meshgrid(lonsraw, latsraw)

        latsraw = []
        lonsraw = []
        tmp.close()

    for files in userVariables.filelist:
        try:
            thisFile = Dataset(files, 'r', format='NETCDF4')
            # clip the dataset according to user lat,lon coordinates
            # mask the data and fill with zeros for later
            tempRaw = thisFile.variables[varName][:, latminIndex:latmaxIndex, lonminIndex:lonmaxIndex].astype('int16')
            tempMask = ma.masked_array(tempRaw, mask=(tempRaw > userVariables.T_BB_MAX), fill_value=0)
            # get the actual values that the mask returned
            # tempMaskedValueOld = ma.zeros((tempRaw.shape)).astype('int16')

            # innerLoopStart = time.time()
            # for index, value in utils.maenumerate(tempMask):
            #    timeIndex, latIndex, lonIndex = index
            #    tempMaskedValueOld[timeIndex, latIndex, lonIndex] = value
            # totalLoopTime+= time.time()-innerLoopStart
            
            # This replaces the loop computation of tempMaskedValueOld above.
            tempMaskedValue = tempMask
            tempMaskedValue[tempMask.mask] = 0
            # Mini unit test of the non-loop version vers loop version. To use,
            # uncomment the loop immediately above
            # assert(np.array_equal(tempMaskedValue,tempMaskedValueOld))

            xtimes = thisFile.variables[timeName]

            # convert this time to a python datastring
            time2store, _ = get_model_times(xtimes, timeName)

            # extend instead of append because get_model_times returns a list already and we don't
            # want a list of list
            timelist.extend(time2store)
            inputData.extend(tempMaskedValue)
            thisFile.close()
            thisFile = None

        except:
            print 'bad file! ', files

    inputData = ma.array(inputData)

    return inputData, timelist, LAT, LON, userVariables
# **********************************************************************************************************************
def get_model_times(xtimes, timeVarName):
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

    baseTimeString = string.lstrip(timeFormat[sinceLoc:])
    baseTime = decode_time_from_string(baseTimeString)

    times = []

    for xtime in xtimes[:]:

        # TODO: This casting may cause problems for data that is hourly with more than one timestep in it
        xtime = int(xtime)

        if xtime % 10 != 0:
            xtime = 1

        if units == 'minutes':
            dt = timedelta(minutes=xtime)
            newTime = baseTime + dt
        elif units == 'hours':
            dt = timedelta(hours=int(xtime))
            newTime = baseTime + dt  # timedelta(hours=int(xtime))
        elif units == 'days':
            dt = timedelta(days=xtime)
            newTime = baseTime + dt
        elif units == 'months':
            # NB. adding months in python is complicated as month length varies and hence ambiguous.+
            # Perform date arithmetic manually
            #  Assumption: the base_date will usually be the first of the month
            #              NB. this method will fail if the base time is on the 29th or higher day of month
            #                      -as can't have, e.g. Feb 31st.
            newMonth = int(baseTime.month + xtime % 12)
            newYear = int(math.floor(baseTime.year + xtime / 12.))
            newTime = datetime.datetime(newYear, newMonth, baseTime.day, baseTime.hour, baseTime.second, 0)
        elif units == 'years':
            dt = datetime.timedelta(years=xtime)
            newTime = baseTime + dt

        times.append(newTime)

    try:
        if len(xtimes) == 1:
            timeStepLength = 0
        else:
            timeStepLength = int(xtimes[1] - xtimes[0] + 1.e-12)

        modelTimeStep = get_model_time_step(units, timeStepLength)

    except:
        raise

    return times, modelTimeStep
# **********************************************************************************************************************
def get_model_time_step(units, stepSize):
    '''
        Purpose:: To determine the time intervals of input data.
                  Levaraged from Apache OCW

        Inputs:: units: a string representing the time units found in the file metadata
                 stepSize: an integer representing the time interval found in the file's metadata

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
        # TODO: need a check for fractional hrs and only one hr i.e. stepSize=0 e.g. with MERG data
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
# **********************************************************************************************************************
def decode_time_from_string(timeString):
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
        if datetime.strptime(timeString, '%Y-%m-%d %H'):
            myTime = datetime.strptime(timeString, '%Y-%m-%d %H')
        elif datetime.strptime(timeString, '%Y-%m-%d %H%M'):
            myTime = datetime.strptime(timeString, '%Y-%m-%d %H%M')
        return myTime

    except ValueError:
        pass

    print 'Error decoding time string: string does not match a predefined time format'
    return 0
    # **********************************************************************************************************************


