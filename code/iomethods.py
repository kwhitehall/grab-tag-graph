import glob
import os
import re
import string
import subprocess
import sys
from datetime import timedelta, datetime

import numpy as np
import numpy.ma as ma
import netCDF4

import utils
import variables

def _get_fileList_for_binaries(dirPath, startTime, endTime):
    '''
        Purpose:: There are two kinds of files in the MERG directory. One has files in binary form, and the other in
                  netCDF. This function is to get only the files that are in between the startTime and endTime that are
                  in binary format.
        Input:: dirPath: Path to where the binary MERG files are.
                startTime: The cutoff for the earliest files you want
                endTime: The cutoff for the latest files you want
        Output:: A list of a path to files that are between startTime and endTime
        Assumptions:: startTime and endTime are in the format YYYYMMDDHHSS
                      The filename always ends with -pixel
                      The date is included in the file name, with the format YYYYMMDDHH
    '''
    fileString = os.path.join(dirPath, '*-pixel')
    fileList = glob.glob(fileString)
    fileList.sort()

    newFileList = []

    start = datetime.strptime(startTime, '%Y%m%d%H%S')  # Convert start and end times to datetimes to compare with
    end = datetime.strptime(endTime, '%Y%m%d%H%S')      # the current file's time

    for file in fileList:
        dateFromFileName = file[5:15]
        dateAsDateTime = datetime.strptime(dateFromFileName, '%Y%m%d%H')

        if start <= end and start <= dateAsDateTime <= end:
            newFileList.append(file)
        elif start <= dateAsDateTime or dateAsDateTime <= end:
            newFileList.append(file)

    return newFileList


def check_for_files(dirPath, startTime, endTime, tdelta, tRes, flag):
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

    if flag is True:
        filelist = filter(os.path.isfile, glob.glob((dirPath+'/*.nc')))
    else:
        filelist = filter(os.path.isfile, glob.glob(dirPath))
    filelist.sort()

    # Check for the filename pattern
    if flag is True:
        for eachPart in re.split(r'[_,-,.,/]', re.split(r'.nc', os.path.basename(filelist[0]))[0]):
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
            fileDate = int(re.search(r'\d+', re.split(r'.nc', os.path.basename(filelist[0]))[0]).group())
            filenamePattern = re.split(str(fileDate), os.path.basename(filelist[0]))[0]
    else:
        for eachPart in re.split(r'[_,-,.,/]', re.split(os.path.basename(filelist[0]))[0]):
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
            fileDate = int(re.search(r'\d+', re.split(os.path.basename(filelist[0]))[0]).group())
            filenamePattern = re.split(str(fileDate), os.path.basename(filelist[0]))[0]

    startFile = glob.glob(dirPath+'/'+filenamePattern + '*' + startTimeInFile)[0]
    endTimeInFile = _find_time_in_file(endTime, startTimeInFile)
    endFile = glob.glob(dirPath+'/'+filenamePattern + '*' + endTimeInFile + '*')[0]

    currFile = startFile
    filelist = []

    # Check for files between startTime and endTime
    while currFile is not endFile:
        if not os.path.isfile(currFile):
            status = False
            return status, filelist
        else:
            filelist.append(currFile)

        status = True
        if currFile == endFile:
            break

        if 'month' in tRes:
            currFileTime += timedelta(days=31*tdelta)
            currTimeInFile = _find_time_in_file(currFileTime.strftime('%Y%m'), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]
        if 'day' in tRes:
            currFileTime += timedelta(days=tdelta)
            currTimeInFile = _find_time_in_file(currFileTime.strftime('%Y%m%d'), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]
        if 'hour' in tRes:
            currFileTime += timedelta(hours=tdelta)
            currTimeInFile = _find_time_in_file(currFileTime.strftime('%Y%m%d%H'), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]
        if 'minute' in tRes:
            currFileTime += timedelta(minutes=tdelta)
            currTimeInFile = _find_time_in_file(currFileTime.strftime('%Y%m%d%H%M'), startTimeInFile)
            currFile = glob.glob(dirPath+'/'+filenamePattern+'*'+currTimeInFile+'*')[0]

    return status, filelist
# **********************************************************************************************************************
def _find_time_in_file(myTime, myTimeInFile):
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

    # Check validity of time
    while utils.valid_date(userVariables.startDateTime) != True:
        print "Invalid time entered for startDateTime!"

    while utils.valid_date(userVariables.endDateTime) != True:
        print "Invalid time entered for endDateTime!"

    # Check if all the files exists in the MERG and TRMM directories entered
    test, _ = check_for_files(userVariables.DIRS['TRMMdirName'], userVariables.startDateTime, userVariables.endDateTime, 3, 'hour', True)
    if test is False:
        print "Error with files in the TRMM directory entered. Please check your files before restarting. "
        return

    test, userVariables.filelist = check_for_files(userVariables.DIRS['CEoriDirName'], userVariables.startDateTime, userVariables.endDateTime, 1, 'hour', True)

    if test is False:
        print "Error with files in the original MERG directory entered. Please check your files before restarting. "
        return

    # Create main directory and file structure for storing intel
    userVariables.DIRS['mainDirStr'] = _create_main_directory(userVariables.DIRS['mainDirStr'])

    return graphVariables
# **********************************************************************************************************************
def _create_main_directory(mainDirStr):
    '''
        Purpose:: To create the main directory for storing information and
                  the subdirectories for storing information
        Input:: mainDir: string representing the directory for where all
                information generated from the program are to be stored
        Returns:: None
        Outputs:: A file structure where data generated from GTG will be stored
        Assumptions:: The user running the program can write at mainDirStr
    '''
    global MAIN_DIRECTORY

    MAIN_DIRECTORY = mainDirStr
    # If directory doesn't exist, create it
    if not os.path.exists(MAIN_DIRECTORY):
        os.makedirs(MAIN_DIRECTORY)

    os.chdir((MAIN_DIRECTORY))
    # Create the subdirectories
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
def read_data(varName, latName, lonName, userVariables, fileType):
    '''
        Purpose::
            Read gridded data into (t, lat, lon) arrays for processing
        Inputs::
            varName: a string representing the variable name to use from the file
            latName: a string representing the latitude from the file's metadata
            lonName: a string representing the longitude from the file's metadata
            userVariables: a UserVariables object
            fileType: a string representing whether we want to read either a netCDF file or a binary file.
        Outputs::
            A 3D masked array (t,lat,lon) with only the variables which meet the minimum temperature
            criteria for each frame if the fileType is netCDF
            If the fileType is binary, then the data is not masked by temperature
        Assumptions::
            (1) All the files requested to extract data are from the same instrument/model, and thus have the same
            metadata properties (varName, latName, lonName) as entered
            (2) Assumes rectilinear grids for input datasets i.e. lat, lon will be 1D arrays
    '''

    # Get list of files between startDateTime and endDateTime in the passed in userVariables object
    if fileType == 'binary':
        userVariables.filelist = _get_fileList_for_binaries(userVariables.DIRS['CEoriDirName'], userVariables.startDateTime,
                                                           userVariables.endDateTime)
    outputData = []
    timelist = []

    # Crash nicely if there are no files
    if len(userVariables.filelist) == 0:
        print 'Error: no files found'
        sys.exit()

    if fileType == 'netCDF':
        tmp = netCDF4.Dataset(userVariables.filelist[0], 'r+', format='NETCDF4')

        alllatsraw = tmp.variables[latName][:]
        alllonsraw = tmp.variables[lonName][:]
        alllonsraw[alllonsraw > 180] = alllonsraw[alllonsraw > 180] - 360.  # convert to -180,180 if necessary

        tmp.close()
    elif fileType == 'binary':
        # Generate latitudes and longitudes because they didn't come with the binary file.
        # See http://www.cpc.ncep.noaa.gov/products/global_precip/html/README for reasoning behind the numbers
        alllatsraw = np.arange(59.982, -60., -0.036383683, dtype=np.float)
        alllonsraw = np.arange(0.0182, 360., 0.036378335, dtype=np.float)

    # Get the lat/lon info data (different resolution)
    latminNETCDF = utils.find_nearest(alllatsraw, float(userVariables.LATMIN))
    latmaxNETCDF = utils.find_nearest(alllatsraw, float(userVariables.LATMAX))
    lonminNETCDF = utils.find_nearest(alllonsraw, float(userVariables.LONMIN))
    lonmaxNETCDF = utils.find_nearest(alllonsraw, float(userVariables.LONMAX))
    latminIndex = (np.where(alllatsraw == latminNETCDF))[0][0]
    latmaxIndex = (np.where(alllatsraw == latmaxNETCDF))[0][0]
    lonminIndex = (np.where(alllonsraw == lonminNETCDF))[0][0]
    lonmaxIndex = (np.where(alllonsraw == lonmaxNETCDF))[0][0]

    # Subset the data
    latsraw = alllatsraw[latminIndex:latmaxIndex]
    lonsraw = alllonsraw[lonminIndex:lonmaxIndex]

    LON, LAT = np.meshgrid(lonsraw, latsraw)

    timeName = 'time'
    for file in userVariables.filelist:
        if fileType == 'netCDF':
            try:
                thisFile = netCDF4.Dataset(file, 'r', format='NETCDF4')

                # Clip the dataset according to user lat, lon coordinates
                # Mask the data and fill with zeros for later
                tempRaw = thisFile.variables[varName][:, latminIndex:latmaxIndex, lonminIndex:lonmaxIndex].astype('int16')
                tempMask = ma.masked_array(tempRaw, mask=(tempRaw > userVariables.T_BB_MAX), fill_value=-999)
                # Get the actual values that the mask returned

                tempMaskedValue = tempMask
                tempMaskedValue[tempMask.mask] = 0

                xtimes = thisFile.variables[timeName]

                # Convert this time to a python datestring
                time2store, _ = _get_model_times(xtimes)

                # Extend instead of append because get_model_times returns a list and we don't want a list of list
                timelist.extend(time2store)
                outputData.extend(tempMaskedValue)
                thisFile.close()
            except:
                print 'bad file!', file

        elif fileType == 'binary':
            try:
                # The data is temperature in Kelvin, for more information
                # go to http://www.cpc.ncep.noaa.gov/products/global_precip/html/README
                data = _read_merg_file(file, shape=(2, 3298, 9896), offset=75.)

                # Subset the data
                data = data[:, latminIndex:latmaxIndex, lonminIndex:lonmaxIndex]

                # Parse time from the file name. Assumes the filename is in the format 'merg_xxxxxxxxxx_4km-pixel'
                dateFromFileName = file[5:15]
                dateAsDateTime = datetime.strptime(dateFromFileName, '%Y%m%d%H')

                timelist.append(dateAsDateTime)
                outputData.extend(data)
            except:
                print 'bad file!', file

    outputData = ma.array(outputData)

    return outputData, timelist, LAT, LON, userVariables
# **********************************************************************************************************************
def _get_model_times(xtimes):
    '''
    Purpose:: Routine to convert from model times ('hours since 1900...', 'days since ...')
    into a python datetime structure. Leveraged from Apache OCW
    Inputs::
        xtimes: a netcdf4.variable object with the times, encoded with the
            described model time format.
    Returns::
        times: a list of python datetime objects describing model data times
        modelTimeStep: a string representing the time step found in the file e.g.
        'hourly','daily','monthly','annual'
    Outputs:: None
    Assumptions:: Assumes the units are described in the variable's "units"
        attribute as a string.
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
    baseTime = _decode_time_from_string(baseTimeString)

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
            # Assumption: the base_date will usually be the first of the month
            #              NB. this method will fail if the base time is on the 29th or higher day of month
            #                      -as can't have, e.g. Feb 31st.
            newMonth = int(baseTime.month + xtime % 12)
            newYear = int(baseTime.year + xtime / 12.)
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

        modelTimeStep = _get_model_time_step(units, timeStepLength)

    except:
        raise

    return times, modelTimeStep
# **********************************************************************************************************************
def _get_model_time_step(units, stepSize):
    '''
        Purpose:: To determine the time intervals of input data.
                  Leveraged from Apache OCW
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
def _decode_time_from_string(timeString):
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
def write_merg_to_ncdf(lonDict, latDict, timeDict, ch4Dict, fileName, dirName, globalAttrDict, dimensionsDict):
    '''
        Purpose:: Write merg data to netCDF format. See the method
                  _read_merg_file to see how the data is arranged.

        Inputs:: lonDict, latDict, timeDict, ch4Dict are all dictionaries with
                     string keys of 'name', 'dataType', and 'dimensions' and their corresponding data as the value.
                     The other key: value pairs are the variable's attributes, except for the last pair which is the
                     data associated with the variable.
                 globalAttrDict is a dictionary where the key is a string that represents the name of the global attribute
                    and the value is the description of the global attribute.
                 dimensionsDict is a dictionary where the key is a string that represents the name of the dimension
                    and the value is either an integer representing the size or none to represent an unlimited dimension
                 fileName is a string representing the fileName where you got the data
                 dirName is a string representing where you want to place the newly created file

        Returns:: newFilePath, which is the path to the new file. It also writes to a netCDF file
        Assumptions:: fileName is in the format: merg_YYYYMMDDHH_4km-pixel

    '''

    newFilePath = os.path.join(dirName, fileName + '.nc')
    ncdf = netCDF4.Dataset(newFilePath, "w", format="NETCDF4")

    # Set global attributes, dimensions, and lastly, each variable
    ncdf.setncatts(globalAttrDict)

    for nameOfDimension, size in dimensionsDict.iteritems():
        ncdf.createDimension(nameOfDimension, size)

    ncdf.createVariable(lonDict['name'], lonDict['dataType'], lonDict['dimensions'])
    ncdf.variables[lonDict['name']].setncattr('units', lonDict['units'])
    ncdf.variables[lonDict['name']].setncattr('long_name', lonDict['long_name'])

    ncdf.createVariable(latDict['name'], latDict['dataType'], latDict['dimensions'])
    ncdf.variables[latDict['name']].setncattr('units', latDict['units'])
    ncdf.variables[latDict['name']].setncattr('long_name', latDict['long_name'])

    ncdf.createVariable(timeDict['name'], timeDict['dataType'], timeDict['dimensions'])

    # Parse time from filename. Assumes the filename is in the format 'merg_xxxxxxxxxx_4km-pixel'
    dateFromFileName = fileName[5:15]
    dateAsDateTime = datetime.strptime(dateFromFileName, '%Y%m%d%H')

    ncdf.variables['time'].setncattr('units', 'hours since ' + str(dateAsDateTime.year) + '-' + str(dateAsDateTime.month) +
                                     '-' + str(dateAsDateTime.day) + ' ' + str(dateAsDateTime.hour))

    ncdf.createVariable(ch4Dict['name'], ch4Dict['dataType'], ch4Dict['dimensions'])
    ncdf.variables[ch4Dict['name']].setncattr('long_name', ch4Dict['long_name'])
    ncdf.variables[ch4Dict['name']].setncattr('time_statistic', ch4Dict['time_statistic'])
    ncdf.variables[ch4Dict['name']].setncattr('missing_value', ch4Dict['missing_value'])

    ncdf.variables['longitude'][:] = lonDict['values']
    ncdf.variables['latitude'][:] = latDict['values']
    ncdf.variables['ch4'][:,:,:] = ch4Dict['values']

    ncdf.close()

    return newFilePath


def _read_merg_file(filepath, shape, offset):
    '''
        Purpose:: Read merg data from a binary file and put it into a numPy array

        Input:: filepath - The path to the merg binary file
                Shape - The shape we want the data to be in
                Offset - Data is sometimes scaled to fit by reducing the value, so add it back if necessary

        Output:: data - A numPy array with datatype np.float with the shape that was passed in

    '''

    mergFile = open(filepath, 'rb')
    mergArray = np.fromfile(mergFile, dtype=np.uint8, count=-1)  # count=-1 means read entire file
    data = mergArray.astype(np.float).reshape(shape)
    data += offset

    mergFile.close()

    return data


def read_netCDF_to_array(filepath, variableToExtract, minT, maxT, minLat, maxLat, minLon, maxLon):
    '''
        Purpose:: Extract the data from a single variable on a netCDF file (from a supported source). The user specifies a
         slice of the data in three dimensions: time, latitude, and longitude. The function returns the data points
         inside this slice. If the slice goes outside of the boundaries of the data available in the file (i.e. if the requested
         3D slice spans a space outside the boundaries of the data file), an exception is thrown.


        Input:: filepath: path to the netCDF file
                filetype: string containing the data source (used to infer the structure, supported: "trmm", "mtsat", wrf")
                variableToExtract: the name of the netCDF variable to be returned
                minT: time of the earliest boundary of the slice request
                maxT: time of the latest boundary of the slice request
                minLat: latitude of the southern boundary of the slice request
                maxLat: latitude of the northern boundary of the slice request
                minLon: longitude of the western boundary of the slice request
                maxLon: longitude of the eastern boundary of the slice request

        Output:: trimmedData: Masked numPy 3D array ([time, lat, lon]) that contains the data requested.
                 trimmedTimes: list of all positions in the time dimension (in 'datetime' object format)
                 trimmedLat: numpy array of all positions in the latitude dimension
                 trimmedLon: numpy array of all positions in the longitude dimension


        Assumption(s):: The variable has only 2 or 3 dimensions. If it has 2, then a third dimension is added.
                        The variable has the dimensions: [time,] latitude, and longitude, in that order

    '''

    # Verifying input
    while -90. > minLat or minLat > 90.:
        try:
            minLat = float(raw_input("minLat passed to read_netCDF_to_array invalid, enter replacement:"))
        except ValueError as e:
            print "Error parsing coordinates: " + e.message
            print "Try again"
    while -90. > maxLat or maxLat > 90.:
        try:
            maxLat = float(raw_input("maxLat passed to read_netCDF_to_array invalid, enter replacement:"))
        except ValueError as e:
            print "Error parsing coordinates: " + e.message
            print "Try again"
    while -180. > minLon or minLon > 180.:
        try:
            minLon = float(raw_input("minLon passed to read_netCDF_to_array invalid, enter replacement:"))
        except ValueError as e:
            print "Error parsing coordinates: " + e.message
            print "Try again"
    while -180. > maxLon or maxLon > 180.:
        try:
            maxLon = float(raw_input("maxLon passed to read_netCDF_to_array invalid, enter replacement:"))
        except ValueError as e:
            print "Error parsing coordinates: " + e.message
            print "Try again"

    dataset = netCDF4.Dataset(filepath, 'r', format='NETCDF4')

    extractedVariable = dataset.variables[variableToExtract][:, :]

    # adding time dimension to 2D data formats
    if extractedVariable.ndim is 2:
        extractedVariable = extractedVariable[np.newaxis, :, :]

    # Finding variables for required dimensions
    variableDimensions = dataset.variables[variableToExtract].dimensions
    if all([(dim in dataset.variables) for dim in variableDimensions]):
        timeVariable = variableDimensions[-3] if len(variableDimensions) >= 3 else None
        latsVariable, lonsVariable = variableDimensions[-2:]
    else:
        timeVariable = None
        lonsVariable, latsVariable = dataset.variables[variableToExtract].coordinates.split()

    if timeVariable == None:
        timeVariable = "Times"

    # Grabbing list of latitudes and longitudes
    latsData = dataset.variables[latsVariable]
    lonsData = dataset.variables[lonsVariable]
    if len(latsData.shape) == 1: # 1D array from data source.
        latsList = latsData[:]
        lonsList = lonsData[:]
    else: # 3D array from data source.
        latsList = latsData[0,:,0]
        lonsList = lonsData[0,0,:]

    # Reformatting longitude format to [-180, 180]
    lonsList[lonsList > 180] = lonsList[lonsList > 180] - 360.

    if dataset.variables is None or timeVariable not in dataset.variables:
        filename = os.path.basename(filepath)

        reObj = re.compile("([0-9]{4})([0-9]{2})([0-9]{2})_([0-9]{2})([0-9]{2})")
        timeMatchGroup = reObj.search(filename)
        timeNumbers = [int(num) for num in timeMatchGroup.groups()]
        time = datetime(*timeNumbers)
        times = [time]
    elif dataset.variables[timeVariable].dtype == np.dtype('S1'):
        stringTimes = ["".join(x) for x in dataset.variables[timeVariable][:]]
        times = [datetime.strptime(x, '%Y-%m-%d_%H:%M:%S') for x in stringTimes]
    elif hasattr(dataset.variables[timeVariable], "units"):
        times = _get_model_times(dataset.variables[timeVariable])[0]
    else:
        print "No approached worked for fetching the time out of the netCDF " \
              "file %s. Enter times." % filepath
        times = []
        for i in range(extractedVariable.shape[0]):
            while True:
                newTimeStr = raw_input("Enter time for data at index %d for "
                                       "variable %s at file %s, in yyyy-mm-dd "
                                       "hh:mm:ss format:"
                                       % (i, variableToExtract, filepath))

                try:
                    newTime = datetime.strptime(newTimeStr, "%Y-%m-%d %H:%M:%S")
                    times.append(newTime)
                    break
                except ValueError as e:
                    print "Error parsing dates: " + e.message
                    print "Try again"

    # Verifying requested area and times are available
    maxLat, maxLon, maxT, minLat, minLon, minT = _check_bounds(latsList,
                                                               lonsList,
                                                               maxLat,
                                                               maxLon,
                                                               maxT,
                                                               minLat,
                                                               minLon,
                                                               minT,
                                                               times)

    # Trimming data according to requested area and times
    trimmedLatsStart = np.where(latsList == utils.find_nearest(latsList, minLat))[0][0]
    trimmedLatsEnd = np.where(latsList == utils.find_nearest(latsList, maxLat))[0][0]

    trimmedLonsStart = np.where(lonsList == utils.find_nearest(lonsList, minLon))[0][0]
    trimmedLonsEnd = np.where(lonsList == utils.find_nearest(lonsList, maxLon))[0][0]

    # Can't use np.where() for times, so use list comprehension instead
    trimmedTimeIndices = [i for i in range(len(times))
                            if times[i] >= minT and times[i] <= maxT]
    trimmedTimeStart = trimmedTimeIndices[0]
    trimmedTimeEnd = trimmedTimeIndices[-1]

    trimmedLats = latsList[trimmedLatsStart:trimmedLatsEnd + 1]
    trimmedLons = lonsList[trimmedLonsStart:trimmedLonsEnd + 1]
    trimmedTimes = times[trimmedTimeStart: trimmedTimeEnd + 1]
    trimmedData = extractedVariable[trimmedTimeStart:trimmedTimeEnd + 1,
                                    trimmedLatsStart:trimmedLatsEnd + 1,
                                    trimmedLonsStart:trimmedLonsEnd + 1]

    dataset.close()

    return ma.masked_array(trimmedData), trimmedTimes, trimmedLats, trimmedLons


def _check_bounds(latsList, lonsList, maxLat, maxLon, maxT, minLat, minLon,
                  minT, times):
    while min(latsList) > minLat or max(latsList) < maxLat:
        print "Requested range is outside file bounds in latitude axis. File " + \
              " bounds are %f, %f. Requested range is %f, %f" % (
                  min(latsList), max(latsList), minLat, maxLat)

        lastMinLat = minLat
        lastMaxLat = maxLat
        try:
            minLat = float(raw_input("Enter new minLat for read_netCDF_to_array:"))
            maxLat = float(raw_input("Enter new maxLat for read_netCDF_to_array:"))
        except ValueError as e:
            print "Error parsing coordinates: " + e.message
            print "Try again"
            minLat = lastMinLat
            maxLat = lastMaxLat

    while min(lonsList) > minLon or max(lonsList) < maxLon:
        print "Requested range is outside file bounds in longitude axis. File " + \
              " bounds are %f, %f. Requested range is %f, %f" % (
                  min(lonsList), max(lonsList), minLon, maxLon)

        lastMinLon = minLon
        lastMaxLon = maxLon

        try:
            minLon = float(raw_input("Enter new minLon for read_netCDF_to_array:"))
            maxLon = float(raw_input("Enter new maxLon for read_netCDF_to_array:"))
        except ValueError as e:
            print "Error parsing coordinates: " + e.message
            print "Try again"
            minLon = lastMinLon
            maxLon = lastMaxLon

    while min(times) > minT or max(times) < maxT:

        print "Requested range is outside file bounds in time axis. File " + \
              " bounds are %s, %s. Requested range is %s, %s" % (
                  min(times), max(times), minT, maxT)
        lastMinT = minT
        lastMaxT = maxT

        minTStr = raw_input("Enter new minT for read_netCDF_to_array in yyyy-mm-dd hh:mm:ss format:")
        maxTStr = raw_input("Enter new maxT for read_netCDF_to_array in yyyy-mm-dd hh:mm:ss format:")

        try:
            minT = datetime.strptime(minTStr, "%Y-%m-%d %H:%M:%S")
            maxT = datetime.strptime(maxTStr, "%Y-%m-%d %H:%M:%S")
        except ValueError as e:
            print "Error parsing dates: " + e.message
            print "Try again"
            minT = lastMinT
            maxT = lastMaxT

    return maxLat, maxLon, maxT, minLat, minLon, minT

# **********************************************************************************************************************


