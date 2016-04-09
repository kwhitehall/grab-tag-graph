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
import netCDF4

import utils
import variables

def get_fileList_for_binaries(dirPath, startTime, endTime):
    '''
        Purpose:: There are two kinds of files in the MERG directory. One has files in binary form, and the other in
                  netCDF. This function is to get only the files that are in between the startTime and endTime that are
                  in binary format.
        Input:: dirPath: Path to where the binary MERG files are.
                startTime: The cutoff for the earliest files you want
                endTime: The cutoff for the latest files you want
        Output:: A list of a path to files that are between
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
        dateFromFileName = [token for token in file.split('_') if token.isdigit()]  # Parse date from MERG binary file name
        dateAsDateTime = datetime.strptime(dateFromFileName[0], '%Y%m%d%H')

        if start <= end and start <= dateAsDateTime <= end:
            newFileList.append(file)
        elif start <= dateAsDateTime or dateAsDateTime <= end:
            newFileList.append(file)

    print len(newFileList)
    return newFileList


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

    # Check for the filename pattern
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
    endFile = glob.glob(dirPath+'/'+filenamePattern + '*' + endTimeInFile + '*')[0]

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

    # Check validity of time
    while utils.valid_date(userVariables.startDateTime) != True:
        print "Invalid time entered for startDateTime!"

    while utils.valid_date(userVariables.endDateTime) != True:
        print "Invalid time entered for endDateTime!"

    # Check if all the files exists in the MERG and TRMM directories entered
    test, _ = check_for_files(userVariables.DIRS['TRMMdirName'], userVariables.startDateTime, userVariables.endDateTime, 3, 'hour')
    if test is False:
        print "Error with files in the TRMM directory entered. Please check your files before restarting. "
        return

    test, userVariables.filelist = check_for_files(userVariables.DIRS['CEoriDirName'], userVariables.startDateTime, userVariables.endDateTime, 1, 'hour')

    if test is False:
        print "Error with files in the original MERG directory entered. Please check your files before restarting. "
        return

    # Create main directory and file structure for storing intel
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

    if fileType == 'binary':
        userVariables.filelist = get_fileList_for_binaries(userVariables.DIRS['CEoriDirName'], userVariables.startDateTime,
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
        alllatsraw, alllonsraw, _ = read_MERG_pixel_file(userVariables.filelist[0])

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
                time2store, _ = get_model_times(xtimes)

                # Extend instead of append because get_model_times returns a list and we don't want a list of list
                timelist.extend(time2store)
                outputData.extend(tempMaskedValue)
                thisFile.close()
            except:
                print 'bad file!', file

        elif fileType == 'binary':
            try:
                _, _, temperatures = read_MERG_pixel_file(file)

                temperatures = temperatures[:, latminIndex:latmaxIndex, lonminIndex:lonmaxIndex]

                dateFromFileName = [token for token in file.split('_') if token.isdigit()]  # Parse date from MERG binary file name
                dateAsDateTime = datetime.strptime(dateFromFileName[0], '%Y%m%d%H')

                timelist.append(dateAsDateTime)
                outputData.extend(temperatures)
            except:
                print 'bad file!', file

    outputData = ma.array(outputData)

    return outputData, timelist, LAT, LON, userVariables
# **********************************************************************************************************************
def get_model_times(xtimes):
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

        modelTimeStep = get_model_time_step(units, timeStepLength)

    except:
        raise

    return times, modelTimeStep
# **********************************************************************************************************************
def get_model_time_step(units, stepSize):
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
        elif datetime.strptime(timeString, '%Y-%m-%d %H%M'):
            myTime = datetime.strptime(timeString, '%Y-%m-%d_%H:%M:%S')
        return myTime

    except ValueError:
        pass

    print 'Error decoding time string: string does not match a predefined time format'
    return 0
    # **********************************************************************************************************************
def write_MERG_pixel_to_ncdf(lonDict, latDict, timeDict, ch4Dict, fileName, dirName, globalAttrDict, dimensionsDict):
    '''
        Purpose:: Write temperature data from specific NumPy arrays to netCDF format. See the method
                  read_MERG_pixel_file to see how the data is arranged.

        Inputs:: lonDict, latDict, timeDict, ch4Dict are all dictionaries with
                     string keys of 'name', 'dataType', and 'dimensions' and their corresponding data as the value.
                     The other key: value pairs are the variable's attributes, except for the last pair which is the
                     data associated with the variable.
                 globalAttrDict is a dictionary where the key is a string that represents the name of the global attribute
                    and the value is the description of the global attribute.
                 dimensionsDict is a dictionary where the key is a string that represents the name of the dimension
                    and the value is either an integer representing the size or none to represent an unlimited dimension
                 fileName is a string representing the fileName where you got the data from
                 dirName is a string representing where you want to place the newly created file

        Returns:: None. It writes to a netCDF file with extension .nc
        Assumptions:: fileName is in the format: merg_YYYYMMDDHH_4km-pixel

    '''

    newFilePath = os.path.join(dirName, fileName + '.nc')
    ncdf = netCDF4.Dataset(newFilePath, "w", format="NETCDF4")

    ncdf.setncatts(globalAttrDict)  # Set global attributes, dimensions, and lastly, each variable

    for nameOfDimension, size in dimensionsDict.iteritems():
        ncdf.createDimension(nameOfDimension, size)

    ncdf.createVariable(lonDict['name'], lonDict['dataType'], lonDict['dimensions'])
    ncdf.variables[lonDict['name']].setncattr('units', lonDict['units'])
    ncdf.variables[lonDict['name']].setncattr('long_name', lonDict['long_name'])

    ncdf.createVariable(latDict['name'], latDict['dataType'], latDict['dimensions'])
    ncdf.variables[latDict['name']].setncattr('units', latDict['units'])
    ncdf.variables[latDict['name']].setncattr('long_name', latDict['long_name'])

    ncdf.createVariable(timeDict['name'], timeDict['dataType'], timeDict['dimensions'])

    dateFromFileName = [token for token in newFilePath.split('_') if token.isdigit()]  # Parse date from MERG binary file name
    dateAsDateTime = datetime.strptime(dateFromFileName[0], '%Y%m%d%H')         # then set attribute for 'time'

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


def read_MERG_pixel_file(path, shape=(2, 3298, 9896), offset=75.):
    '''
        Purpose:: Read MERG brightness temperature from binary file. Thanks to Brian Wilson for this contribution.
                  File contains two large arrays (2 time epochs: on the hour and the half hour)
                  of temperature (Kelvin) as an unsigned integer byte, offset by 75 so it will fit in the 0-255 range.

                  For documentation, see http://www.cpc.ncep.noaa.gov/products/global_precip/html/README

        Input:: Path - The path to the MERG binary file
                Shape - The shape we want the data to be in
                Offset - The temperatures were scaled to fit into 1-byte by subtracting 75, so the offset is used to
                         add 75 back to each temperature.
        Output:: lon - A numPy array containing longitudes from 0 to 360 degrees
                 lat - A numpy array containing latitudes from 60 to -60 degrees

        Assumption::The binary file was unmodified when downloaded. The shape tuple that is hardcoded in the parameter
                    will always be the same unless the data changes.

    '''
    pixel_file = open(path, 'rb')
    pixel_array = np.fromfile(pixel_file, dtype=np.uint8, count=-1)  # count=-1 means read entire file
    temperatures = pixel_array.astype(np.float).reshape(shape)
    temperatures += offset

    lon = np.arange(0.0182, 360., 0.036378335, dtype=np.float)
    lat = np.arange(59.982, -60., -0.036383683, dtype=np.float)

    pixel_file.close()

    return lon, lat, temperatures


# TODO this seems to be unnecessary but I'll leave it for now just in case
FORMAT_DEFS = {
    "trmm" :
        {
            "longitude" : "longitude",
            "latitude" : "latitude",

            "time_handling": {
                "method": "get_model_times",
                "variable": "time"
            }
        },
    "mtsat" :
        {
            "longitude": "longitude",
            "latitude": "latitude",

            "time_handling": {
                "method" :"filename_time_regex",
                "regex_object": re.compile("([0-9]{4})([0-9]{2})([0-9]{2})_([0-9]{2})([0-9]{2})")
            }
        },
    "wrf":
        {
            "longitude": "XLONG",
            "latitude": "XLAT",

            "time_handling": {
                "method": "string_times"
            }
        }
}

def read_netCDF_to_array(filepath, filetype, variable_to_extract, min_lat, max_lat, min_lon, max_lon,
                         min_t, max_t):
    '''
        TODO add documentation
    '''

    dataset = netCDF4.Dataset(filepath, 'r', format='NETCDF4')

    extracted_variable = dataset.variables[variable_to_extract][:, :]

    # adding time dimension to 2D data formats
    if extracted_variable.ndim is 2:
        extracted_variable = extracted_variable[np.newaxis, :, :]

    lats_data = dataset.variables[FORMAT_DEFS[filetype]["latitude"]]
    lons_data = dataset.variables[FORMAT_DEFS[filetype]["longitude"]]

    # Grabing list of latitudes and longitudes
    if filetype is "wrf":
        # TODO slicing assumes dimentions are always in this order : u'Time', u'south_north', u'west_east'
        # TODO check if this is true
        lats_list = lats_data[0,:,0]
        lons_list = lons_data[0,0,:]
    else:
        lats_list = lats_data[:]
        lons_list = lons_data[:]

    # Reformatting longitude format to [-180, 180]
    lons_list[lons_list > 180] = lons_list[lons_list > 180] - 360.

    # Grabbing list of times
    time_dict = FORMAT_DEFS[filetype]["time_handling"]
    if time_dict["method"] is "filename_time_regex":
        filename = path.basename(filepath)

        time_match_group = time_dict["regex_object"].search(filename)
        time_numbers = [int(num) for num in time_match_group.groups()]
        time = datetime(*time_numbers)
        times = [time]
    elif time_dict["method"] is "get_model_times":
        times = get_model_times(dataset.variables[time_dict["variable"]])
    elif time_dict["method"] is "string_times":
        string_times = ["".join(x) for x in dataset.variables["Times"][:]]
        times = [decode_time_from_string(x) for x in string_times]

    # Verifying requested area and times are available
    if min(lats_list) > min_lat or min(lons_list) > min_lon or \
            max(lats_list) < max_lat or max(lons_list) < max_lon:
        raise RuntimeError("Area requested is outside of file data range")
    if min(times) > min_t or max(times) < max_t:
        raise RuntimeError("Time slice requested is outside of file data range")

    # Trimming data according to requested area and times
    trimmed_lats_indices = [i for i in range(len(lats_list))
                        if lats_list[i] >= min_lat and lats_list[i] <= max_lat]
    trimmed_lats_start = trimmed_lats_indices[0]
    trimmed_lats_end = trimmed_lats_indices[-1] + 1

    trimmed_lons_indices = [i for i in range(len(lons_list))
                        if lons_list[i] >= min_lon and lons_list[i] <= max_lon]
    trimmed_lons_start = trimmed_lons_indices[0]
    trimmed_lons_end = trimmed_lons_indices[-1] + 1

    trimmed_time_indices = [i for i in range(len(times))
                        if times >= min_t and times <= max_t]
    trimmed_time_start = trimmed_time_indices[0]
    trimmed_time_end = trimmed_time_indices[-1] + 1


    trimmed_lats = lats_list[trimmed_lats_start:trimmed_lats_end]
    trimmed_lons = lons_list[trimmed_lons_start:trimmed_lons_end]
    trimmed_times = times[trimmed_time_start: trimmed_time_end]

    trimmed_data = extracted_variable[trimmed_time_start: trimmed_time_end,
                                      trimmed_lats_start:trimmed_lats_end,
                                      trimmed_lons_start:trimmed_lons_end]


    dataset.close()

    return ma.masked_array(trimmed_data), trimmed_times, trimmed_lats, trimmed_lons


# **********************************************************************************************************************



