from datetime import timedelta, datetime
import itertools
import os
import subprocess
import re
import string
import math

from netCDF4 import Dataset
import numpy.ma as ma
import numpy as np
from scipy.ndimage import map_coordinates


def checkForFiles(startTime, endTime, thisDir, fileType):
    '''
    Purpose:: To ensure all the files between the starttime and endTime
              exist in the directory supplied

    Input:: 
            startTime: a string yyyymmmddhh representing the starttime 
            endTime: a string yyyymmmddhh representing the endTime
            thisDir: a string representing the directory path where to 
                look for the file
            fileType: an integer representing the type of file in the directory
                1 - MERG original files, 2 - TRMM original files

    Output:: 
            status: a boolean representing whether all files exists

    '''
    filelist=[]
    startFilename = ''
    endFilename =''
    currFilename = ''
    status = False
    startyr = int(startTime[:4])
    startmm = int(startTime[4:6])
    startdd = int(startTime[6:8])
    starthr = int(startTime[-2:])
    endyr = int(endTime[:4])
    endmm = int(endTime[4:6])
    enddd = int(endTime[6:8])
    endhh = int(endTime[-2:])
    curryr = startyr
    currmm = startmm
    currdd = startdd
    currhr = starthr
    currmmStr = ''
    currddStr = ''
    currhrStr = ''
    endmmStr = ''
    endddStr =''
    endhhStr = ''

    #check that the startTime is before the endTime
    if fileType == 1:
        startFilename = "merg_"+startTime+"_4km-pixel.nc"
        endFilename = thisDir+"/merg_"+endTime+"_4km-pixel.nc"

    if fileType == 2:
        #TODO:: determine closest time for TRMM files for end 
        #http://disc.sci.gsfc.nasa.gov/additional/faq/precipitation_faq.shtml#convert
        if starthr%3 == 2:
            currhr += 1 
        elif starthr%3 ==1:
            currhr -= 1
        else:
            currhr = starthr

        curryr, currmmStr, currddStr, currhrStr,_,_,_ = findTime(curryr, currmm, currdd, currhr)

        startFilename = "3B42."+str(curryr)+currmmStr+currddStr+"."+currhrStr+".7A.nc"  
        if endhh%3 == 2:
            endhh += 1
        elif endhh%3 ==1:
            endhh -= 1

        endyr, endmmStr, endddStr, endhhStr, _, _, _ = findTime(endyr, endmm, enddd, endhh)

        endFilename = thisDir+"/3B42."+str(endyr)+endmmStr+endddStr+"."+endhhStr+".7A.nc"

    #check for files between startTime and endTime
    currFilename = thisDir+"/"+startFilename

    while currFilename is not endFilename:
        if not os.path.isfile(currFilename):
            print "file is missing! Filename: ", currFilename
            status = False
            return status,filelist
        else:
            #create filelist
            filelist.append(currFilename)
        
        status = True
        if currFilename == endFilename:
            break

        #generate new currFilename
        if fileType == 1:
            currhr +=1
        elif fileType ==2:
            currhr += 3

        curryr, currmmStr, currddStr, currhrStr, currmm, currdd, currhr = findTime(curryr, currmm, currdd, currhr)

        if fileType == 1:
            currFilename = thisDir+"/"+"merg_"+str(curryr)+currmmStr+currddStr+currhrStr+"_4km-pixel.nc"
        if fileType == 2:
            currFilename = thisDir+"/"+"3B42."+str(curryr)+currmmStr+currddStr+"."+currhrStr+".7A.nc"

    return status,filelist
#******************************************************************
def createMainDirectory(mainDirStr):
    '''
    Purpose:: 
        To create the main directory for storing information and
        the subdirectories for storing information
    Input:: 
        mainDir: a directory for where all information generated from
            the program are to be stored
    Output:: None

    '''
    global MAINDIRECTORY

    MAINDIRECTORY = mainDirStr
    #if directory doesnt exist, creat it
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
def decodeTimeFromString(time_string):
    '''
    Taken from process.py
     Decodes string into a python datetime object
     *Method:* tries a bunch of different time format possibilities and hopefully one of them will hit.
     ::

       **Input:**  time_string - a string that represents a date/time

       **Output:** mytime - a python datetime object
    '''
    # This will deal with times that use decimal seconds
    if '.' in time_string:
        time_string = time_string.split('.')[0] + '0'
    else:
        pass

    try:
        mytime = datetime.strptime(time_string,'%Y-%m-%d %H')
        return mytime

    except ValueError:
        pass

    print 'Error decoding time string: string does not match a predefined time format'
    return 0
#******************************************************************
def doRegrid(q, lat, lon, lat2, lon2, order=1, mdi= -999999999):
    '''
     Perform regridding from one set of lat,lon values onto a new set (lat2,lon2)
    
     Input::
         q          - the variable to be regridded
         lat,lon    - original co-ordinates corresponding to q values
         lat2,lon2  - new set of latitudes and longitudes that you want to regrid q onto 
         order      - (optional) interpolation order 1=bi-linear, 3=cubic spline
         mdi        - (optional) fill value for missing data (used in creation of masked array)
      
     Output::
         q2  - q regridded onto the new set of lat2,lon2 
    
    '''

    nlat = q.shape[0]
    nlon = q.shape[1]

    nlat2 = lat2.shape[0]
    nlon2 = lon2.shape[1]

    # To make our lives easier down the road, let's 
    # turn these into arrays of x & y coords
    loni = lon2.ravel()
    lati = lat2.ravel()

    loni = loni.copy() # NB. it won't run unless you do this...
    lati = lati.copy()

    # Now, we'll set points outside the boundaries to lie along an edge
    loni[loni > lon.max()] = lon.max()
    loni[loni < lon.min()] = lon.min()
    
    # To deal with the "hard" break, we'll have to treat y differently,
    # so we're just setting the min here...
    lati[lati > lat.max()] = lat.max()
    lati[lati < lat.min()] = lat.min()
    
    
    # We need to convert these to (float) indicies
    #   (xi should range from 0 to (nx - 1), etc)
    loni = (nlon - 1) * (loni - lon.min()) / (lon.max() - lon.min())
    
    # Deal with the "hard" break in the y-direction
    lati = (nlat - 1) * (lati - lat.min()) / (lat.max() - lat.min())
    
    # Notes on dealing with MDI when regridding data.
    #  Method adopted here:
    #    Use bilinear interpolation of data by default (but user can specify other order using order=... in call)
    #    Perform bilinear interpolation of data, and of mask.
    #    To be conservative, new grid point which contained some missing data on the old grid is set to missing data.
    #            -this is achieved by looking for any non-zero interpolated mask values.
    #    To avoid issues with bilinear interpolation producing strong gradients leading into the MDI,
    #     set values at MDI points to mean data value so little gradient visible = not ideal, but acceptable for now.
    
    # Set values in MDI so that similar to surroundings so don't produce large gradients when interpolating
    # Preserve MDI mask, by only changing data part of masked array object.
    for shift in (-1, 1):
        for axis in (0, 1):
            q_shifted = np.roll(q, shift=shift, axis=axis)
            idx = ~q_shifted.mask * q.mask
            q.data[idx] = q_shifted[idx]

    # Now we actually interpolate
    # map_coordinates does cubic interpolation by default, 
    # use "order=1" to preform bilinear interpolation instead...
    q2 = map_coordinates(q, [lati, loni], order=order)
    q2 = q2.reshape([nlat2, nlon2])

    # Set values to missing data outside of original domain
    q2 = ma.masked_array(q2, mask=np.logical_or(np.logical_or(lat2 >= lat.max(), 
                                                              lat2 <= lat.min()), 
                                                np.logical_or(lon2 <= lon.min(), 
                                                              lon2 >= lon.max())))
    
    # Make second map using nearest neighbour interpolation -use this to determine locations with MDI and mask these
    qmdi = np.zeros_like(q)
    qmdi[q.mask == True] = 1.
    qmdi[q.mask == False] = 0.
    qmdi_r = map_coordinates(qmdi, [lati, loni], order=order)
    qmdi_r = qmdi_r.reshape([nlat2, nlon2])
    mdimask = (qmdi_r != 0.0)
    
    # Combine missing data mask, with outside domain mask define above.
    q2.mask = np.logical_or(mdimask, q2.mask)

    return q2
#******************************************************************
def findNearest(thisArray,value):
    '''
    Purpose :: to determine the value within an array closes to 
            another value

    Input ::
    Output::
    '''
    idx = (np.abs(thisArray-value)).argmin()
    return thisArray[idx]
#****************************************************************** 
def findTime(curryr, currmm, currdd, currhr):
    '''
    Purpose:: To determine the new yr, mm, dd, hr

    Input:: curryr, an integer representing the year
            currmm, an integer representing the month
            currdd, an integer representing the day
            currhr, an integer representing the hour

    Output::curryr, an integer representing the year
            currmm, an integer representing the month
            currdd, an integer representing the day
            currhr, an integer representing the hour
    '''
    if currhr > 23:
        currhr = 0
        currdd += 1
        if currdd > 30 and (currmm == 4 or currmm == 6 or currmm == 9 or currmm == 11):
            currmm +=1
            currdd = 1
        elif currdd > 31 and (currmm == 1 or currmm ==3 or currmm == 5 or currmm == 7 or currmm == 8 or currmm == 10):
            currmm +=1
            currdd = 1
        elif currdd > 31 and currmm == 12:
            currmm = 1
            currdd = 1
            curryr += 1
        elif currdd > 28 and currmm == 2 and (curryr%4)!=0:
            currmm = 3
            currdd = 1
        elif (curryr%4)==0 and currmm == 2 and currdd>29:
            currmm = 3
            currdd = 1

    if currmm < 10:
        currmmStr="0"+str(currmm)
    else:
        currmmStr = str(currmm)

    if currdd < 10:
        currddStr = "0"+str(currdd)
    else:
        currddStr = str(currdd)

    if currhr < 10:
        currhrStr = "0"+str(currhr)
    else:
        currhrStr = str(currhr)

    return curryr, currmmStr, currddStr, currhrStr, currmm, currdd, currhr
#****************************************************************** 
def getModelTimes(xtimes, timeVarName):
    '''
    Taken from process.py, removed the file opening at the beginning 
    TODO:  Do a better job handling dates here
    Routine to convert from model times ('hours since 1900...', 'days since ...')
    into a python datetime structure

    Input::
        modelFile - path to the model tile you want to extract the times list and modelTimeStep from
        timeVarName - name of the time variable in the model file

    Output::
        times  - list of python datetime objects describing model data times
        modelTimeStep - 'hourly','daily','monthly','annual'
    '''

    timeFormat = xtimes.units
    # search to check if 'since' appears in units
    try:
        sinceLoc = re.search('since', timeFormat).end()

    except AttributeError:
        print 'Error decoding model times: time variable attributes do not contain "since"'
        raise

    units = None
    TIME_UNITS = ('minutes', 'hours', 'days', 'months', 'years')
    # search for 'seconds','minutes','hours', 'days', 'months', 'years' so know units
    for unit in TIME_UNITS:
        if re.search(unit, timeFormat):
            units = unit
            break

    # cut out base time (the bit following 'since')
    base_time_string = string.lstrip(timeFormat[sinceLoc:])
    # decode base time
    base_time = decodeTimeFromString(base_time_string)
    
    times = []

    for xtime in xtimes[:]:
        
        if xtime%10 != 0:
            xtime = 1

        # Cast time as an int
        #TODO: KDW this may cause problems for data that is hourly with more than one timestep in it
        xtime = int(xtime) 
        
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
            # Perform date arithmatic manually
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
       
        #if timeStepLength is zero do not normalize times as this would create an empty list for MERG (hourly) data
        if timeStepLength != 0:
            times = normalizeDatetimes(times, modelTimeStep) 
    except:
        raise

    return times, modelTimeStep
#******************************************************************
def getModelTimeStep(units, stepSize):
    # Time units are now determined. Determine the time intervals of input data (mdlTimeStep)
    '''
    Taken from process.py
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
        #need a check for fractional hrs and only one hr i.e. stepSize=0 e.g. with MERG data
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
def maenumerate(mArray):
    '''
    Purpose::
        Utility script for returning the actual values from the masked array
        Taken from: http://stackoverflow.com/questions/8620798/numpy-ndenumerate-for-masked-arrays
    
    Input::
        mArray: the masked array returned from the ma.array() command
        
        
    Output::
        maskedValues: 3D (t,lat,lon), value of only masked values
    
    '''

    mask = ~mArray.mask.ravel()
    #beware yield fast, but generates a type called "generate" that does not allow for array methods
    for index, maskedValue in itertools.izip(np.ndenumerate(mArray), mask):
        if maskedValue: 
            yield index 
#******************************************************************
def preprocessingMERG(MERGdirname):
    '''
    Purpose::
        Utility script for unzipping and converting the merg*.Z files from Mirador to 
        NETCDF format. The files end up in a folder called mergNETCDF in the directory
        where the raw MERG data is
        NOTE: VERY RAW AND DIRTY 

    Input::
        Directory to the location of the raw MERG files, preferably zipped
        
    Output::
       none

    Assumptions::
       1 GrADS (http://www.iges.org/grads/gadoc/) and lats4D (http://opengrads.org/doc/scripts/lats4d/)
         have been installed on the system and the user can access 
       2 User can write files in location where script is being called
       3 the files havent been unzipped 
    '''

    os.chdir((MERGdirname+'/'))
    imgFilename = ''

    #Just incase the X11 server is giving problems
    subprocess.call('export DISPLAY=:0.0', shell=True)

    for files in glob.glob("*-pixel"):
    #for files in glob.glob("*.Z"):
        fname = os.path.splitext(files)[0]

        #unzip it
        #bash_cmd = 'gunzip ' + files
        #subprocess.call(bash_cmd, shell=True)

        #determine the time from the filename
        ftime = re.search('\_(.*)\_',fname).group(1)

        yy = ftime[0:4]
        mm = ftime[4:6]
        day = ftime[6:8]
        hr = ftime [8:10]

        #TODO: must be something more efficient!

        if mm=='01':
            mth = 'Jan'
        if mm == '02':
            mth = 'Feb'
        if mm == '03':
            mth = 'Mar'
        if mm == '04':
            mth = 'Apr'
        if mm == '05':
            mth = 'May'
        if mm == '06':
            mth = 'Jun'
        if mm == '07':
            mth = 'Jul'
        if mm == '08':
            mth = 'Aug'
        if mm == '09':
            mth = 'Sep'
        if mm == '10':
            mth = 'Oct'
        if mm == '11':
            mth = 'Nov'
        if mm == '12':
            mth = 'Dec'


        subprocess.call('rm merg.ctl', shell=True)
        subprocess.call('touch merg.ctl', shell=True)
        replaceExpDset = 'echo DSET ' + fname +' >> merg.ctl'
        replaceExpTdef = 'echo TDEF 99999 LINEAR '+hr+'z'+day+mth+yy +' 30mn' +' >> merg.ctl'
        subprocess.call(replaceExpDset, shell=True) 
        subprocess.call('echo "OPTIONS yrev little_endian template" >> merg.ctl', shell=True)
        subprocess.call('echo "UNDEF  330" >> merg.ctl', shell=True)
        subprocess.call('echo "TITLE  globally merged IR data" >> merg.ctl', shell=True)
        subprocess.call('echo "XDEF 9896 LINEAR   0.0182 0.036378335" >> merg.ctl', shell=True)
        subprocess.call('echo "YDEF 3298 LINEAR   -59.982 0.036383683" >> merg.ctl', shell=True)
        subprocess.call('echo "ZDEF   01 LEVELS 1" >> merg.ctl', shell=True)
        subprocess.call(replaceExpTdef, shell=True)
        subprocess.call('echo "VARS 1" >> merg.ctl', shell=True)
        subprocess.call('echo "ch4  1  -1,40,1,-1 IR BT  (add  "75" to this value)" >> merg.ctl', shell=True)
        subprocess.call('echo "ENDVARS" >> merg.ctl', shell=True)

        #generate the lats4D command for GrADS
        #lats4D = 'lats4d -v -q -lat '+LATMIN + ' ' +LATMAX +' -lon ' +LONMIN +' ' +LONMAX +' -time '+hr+'Z'+day+mth+yy + ' -func @+75 ' + '-i merg.ctl' + ' -o ' + fname
        
        #lats4D = 'lats4d -v -q -lat -40 -15 -lon 10 40 -time '+hr+'Z'+day+mth+yy + ' -func @+75 ' + '-i merg.ctl' + ' -o ' + fname
        #lats4D = 'lats4d -v -q -lat -5 40 -lon -90 60 -func @+75 ' + '-i merg.ctl' + ' -o ' + fname

        #gradscmd = 'grads -blc ' + '\'' +lats4D + '\''
        #run grads and lats4d command
        #subprocess.call(gradscmd, shell=True)
        imgFilename = hr+'Z'+day+mth+yy+'.gif'
        tempMaskedImages(imgFilename)

    #when all the files have benn converted, mv the netcdf files
    #subprocess.call('mkdir mergNETCDF', shell=True)
    #subprocess.call('mv *.nc mergNETCDF', shell=True)
    #mv all images
    subprocess.call('mkdir mergImgs', shell=True)
    subprocess.call('mv *.gif mergImgs', shell=True)
    return
#******************************************************************
def postProcessingNetCDF(dataset, dirName = None):
    '''
    
    TODO: UPDATE TO PICK UP LIMITS FROM FILE FOR THE GRADS SCRIPTS

    Purpose::
        Utility script displaying the data in generated NETCDF4 files 
        in GrADS
        NOTE: VERY RAW AND DIRTY 

    Input::
        dataset: integer representing post-processed MERG (1) or TRMM data (2) or original MERG(3)
        string: Directory to the location of the raw (MERG) files, preferably zipped
        
    Output::
       images in location as specfied in the code

    Assumptions::
       1 GrADS (http://www.iges.org/grads/gadoc/) and lats4D (http://opengrads.org/doc/scripts/lats4d/)
         have been installed on the system and the user can access 
       2 User can write files in location where script is being called  
    ''' 
    
    coreDir = os.path.dirname(os.path.abspath(__file__))
    ImgFilename = ''
    frameList=[]
    fileList =[]
    lines =[]
    var =''
    firstTime = True
    printLine = 0
    lineNum = 1
    #Just incase the X11 server is giving problems
    subprocess.call('export DISPLAY=:0.0', shell=True)

    prevFrameNum = 0

    if dataset == 1:
        var = 'ch4'
        ctlTitle = 'TITLE MCC search Output Grid: Time  lat lon'
        ctlLine = 'brightnesstemp=\>ch4     1  t,y,x    brightnesstemperature'
        origsFile = coreDir+"/../GrADSscripts/cs1.gs"
        gsFile = coreDir+"/../GrADSscripts/cs2.gs"
        sologsFile = coreDir+"/../GrADSscripts/mergeCE.gs"
        lineNum = 50
    
    elif dataset ==2:
        var = 'precipAcc'
        ctlTitle ='TITLE  TRMM MCS accumulated precipitation search Output Grid: Time  lat lon '
        ctlLine = 'precipitation_Accumulation=\>precipAcc     1  t,y,x    precipAccu'
        origsFile = coreDir+"/../GrADSscripts/cs3.gs"
        gsFile = coreDir+"/../GrADSscripts/cs4.gs"
        sologsFile = coreDir+"/../GrADSscripts/TRMMCE.gs"
        lineNum = 10

    elif dataset ==3:
        var = 'ch4'
        ctlTitle = 'TITLE MERG DATA'
        ctlLine = 'ch4=\>ch4     1  t,y,x    brightnesstemperature'
        origsFile = coreDir+"/../GrADSscripts/cs1.gs"
        sologsFile = coreDir+"/../GrADSscripts/infrared.gs"
        lineNum = 54            

    #sort files
    os.chdir((dirName+'/'))
    try:
        os.makedirs('ctlFiles')
    except:
        print "ctl file folder created already"
        
    files = filter(os.path.isfile, glob.glob("*.nc"))
    files.sort(key=lambda x: os.path.getmtime(x))
    for eachfile in files:
        fullFname = os.path.splitext(eachfile)[0]
        fnameNoExtension = fullFname.split('.nc')[0]
        
        if dataset == 2 and fnameNoExtension[:4] != "TRMM":
            continue

        if dataset == 1 or dataset == 2:
            frameNum = int((fnameNoExtension.split('CE')[0]).split('00F')[1])
        
        #create the ctlFile
        ctlFile1 = dirName+'/ctlFiles/'+fnameNoExtension + '.ctl'
        #the ctl file
        subprocessCall = 'rm ' +ctlFile1
        subprocess.call(subprocessCall, shell=True)
        subprocessCall = 'touch '+ctlFile1
        subprocess.call(subprocessCall, shell=True)
        lineToWrite = 'echo DSET ' + dirName+'/'+fnameNoExtension+'.nc' +' >>' + ctlFile1 
        subprocess.call(lineToWrite, shell=True)  
        lineToWrite = 'echo DTYPE netcdf >> '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite = 'echo UNDEF 0 >> '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite = 'echo '+ctlTitle+' >> '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        fname = dirName+'/'+fnameNoExtension+'.nc'
        if os.path.isfile(fname):   
            #open NetCDF file add info to the accu 
            print "opening file ", fname
            fileData = Dataset(fname,'r',format='NETCDF4')
            lats = fileData.variables['latitude'][:]
            lons = fileData.variables['longitude'][:]
            LONDATA, LATDATA = np.meshgrid(lons,lats)
            nygrd = len(LATDATA[:,0]) 
            nxgrd = len(LONDATA[0,:])
            fileData.close()
        lineToWrite = 'echo XDEF '+ str(nxgrd) + ' LINEAR ' + str(min(lons)) +' '+ str((max(lons)-min(lons))/nxgrd) +' >> ' +ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite = 'echo YDEF '+ str(nygrd)+' LINEAR  ' + str(min(lats)) + ' ' + str((max(lats)-min(lats))/nygrd) +' >> '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite = 'echo ZDEF   01 LEVELS 1 >> '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite = 'echo TDEF 99999 linear 31aug2009 1hr >> '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite = 'echo VARS 1 >> '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite ='echo '+ctlLine+' >> '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite = 'echo ENDVARS >>  '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite = 'echo  >>  '+ctlFile1
        subprocess.call(lineToWrite, shell=True)

        #create plot of just that data
        subprocessCall = 'cp '+ origsFile+' '+sologsFile
        subprocess.call(subprocessCall, shell=True)

        ImgFilename = fnameNoExtension + '.gif'
                    
        displayCmd = '\''+'d '+ var+'\''+'\n'
        newFileCmd = '\''+'open '+ ctlFile1+'\''+'\n'
        colorbarCmd = '\''+'run cbarn'+'\''+'\n'
        printimCmd = '\''+'printim '+MAINDIRECTORY+'/images/'+ImgFilename+' x800 y600 white\''+'\n'
        quitCmd = '\''+'quit'+'\''+'\n'
            
        GrADSscript = open(sologsFile,'r+')
        lines1 = GrADSscript.readlines()
        GrADSscript.seek(0)
        lines1.insert((1),newFileCmd)
        lines1.insert((lineNum+1),displayCmd)
        lines1.insert((lineNum+2), colorbarCmd)
        lines1.insert((lineNum+3), printimCmd)
        lines1.insert((lineNum + 4), quitCmd)
        GrADSscript.writelines(lines1)
        GrADSscript.close()
        #run the script
        runGrads = 'run '+ sologsFile
        gradscmd = 'grads -blc ' + '\'' +runGrads + '\''+'\n'
        subprocess.call(gradscmd, shell=True)

        if dataset == 1 or dataset == 2:

            if prevFrameNum != frameNum and firstTime == False:
                #counter for number of files (and for appending info to lines)
                count = 0
                subprocessCall = 'cp '+ origsFile+ ' '+gsFile
                subprocess.call(subprocessCall, shell=True)
                for fileName in frameList:
                    if count == 0:
                        frame1 = int((fileName.split('.nc')[0].split('CE')[0]).split('00F')[1])

                    fnameNoExtension = fileName.split('.nc')[0]
                    frameNum = int((fnameNoExtension.split('CE')[0]).split('00F')[1])
                    
                    if frameNum == frame1: 
                        CE_num = fnameNoExtension.split('CE')[1]
                        ImgFilename = fnameNoExtension.split('CE')[0] + '.gif'
                        ctlFile1 = dirName+'/ctlFiles/'+fnameNoExtension + '.ctl'

                        #build cs.gs script will all the CE ctl files and the appropriate display command
                        newVar = var+'.'+CE_num
                        newDisplayCmd = '\''+'d '+ newVar+'\''+'\n'
                        newFileCmd = '\''+'open '+ ctlFile1+'\''+'\n'
                        GrADSscript = open(gsFile,'r+')
                        lines1 = GrADSscript.readlines()
                        GrADSscript.seek(0)
                        lines1.insert((1+count),newFileCmd)
                        lines1.insert((lineNum+count+1),newDisplayCmd)
                        GrADSscript.writelines(lines1)
                        GrADSscript.close()
                    count +=1

                colorbarCmd = '\''+'run cbarn'+'\''+'\n'
                printimCmd = '\''+'printim '+MAINDIRECTORY+'/images/'+ImgFilename+' x800 y600 white\''+'\n'
                quitCmd = '\''+'quit'+'\''+'\n'
                GrADSscript = open(gsFile,'r+')
                lines1 = GrADSscript.readlines()
                GrADSscript.seek(0)
                lines1.insert((lineNum+(count*2)+1), colorbarCmd)
                lines1.insert((lineNum + (count*2)+2), printimCmd)
                lines1.insert((lineNum + (count*2)+3), quitCmd)
                GrADSscript.writelines(lines1)
                GrADSscript.close()
                
                #run the script
                runGrads = 'run '+ gsFile
                gradscmd = 'grads -blc ' + '\'' +runGrads + '\''+'\n'
                subprocess.call(gradscmd, shell=True)
                
                #remove the file data stuff
                subprocessCall = 'cd '+dirName
                
                #reset the list for the next frame
                fileList = frameList
                frameList = []
                for thisFile in fileList:
                    if int(((thisFile.split('.nc')[0]).split('CE')[0]).split('00F')[1]) == frameNum:
                        frameList.append(thisFile)
                frameList.append(eachfile)
                prevFrameNum = frameNum
                
            else:
                frameList.append(eachfile)
                prevFrameNum = frameNum
                firstTime = False
                
    return  
#****************************************************************** 
def tempMaskedImages(imgFilename):
    '''
    Purpose:: 
        To generate temperature-masked images for a first pass verification

    Input::
        imgFilename: filename for the gif file
        
    Output::
        None - Gif images for each file of T_bb less than 250K are generated in folder called mergImgs

    Assumptions::
       Same as for preprocessingMERG
       1 GrADS (http://www.iges.org/grads/gadoc/) and lats4D (http://opengrads.org/doc/scripts/lats4d/)
         have been installed on the system and the user can access 
       2 User can write files in location where script is being called
       3 the files havent been unzipped 
    '''
    
    subprocess.call('rm tempMaskedImages.gs', shell=True)
    subprocess.call('touch tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'open merg.ctl''\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'set mpdset hires''\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'set lat -5 30''\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'set lon -40 30''\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'set cint 10''\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'set clevs 190 200 210 220 230 240 250''\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'set gxout shaded''\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'d ch4+75''\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'run cbarn''\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'draw title Masked Temp @ '+imgFilename +'\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'printim '+imgFilename +' x1000 y800''\'" >> tempMaskedImages.gs', shell=True)
    subprocess.call('echo "''\'quit''\'" >> tempMaskedImages.gs', shell=True)
    gradscmd = 'grads -blc ' + '\'run tempMaskedImages.gs''\'' 
    subprocess.call(gradscmd, shell=True)
    return
#******************************************************************
def validDate(dataString):
    '''
    '''

    if len(dataString) > 10:
        print "invalid time entered"
        return 0

    yr = int(dataString[:4])
    mm = int(dataString[4:6])
    dd = int(dataString[6:8])
    hh = int(dataString[-2:])

    if mm < 1 or mm > 12:
        return 0
    elif hh < 0 or hh > 23:
        return 0
    elif (dd< 0 or dd > 30) and (mm == 4 or mm == 6 or mm == 9 or mm == 11):
        return 0
    elif (dd< 0 or dd > 31) and (mm == 1 or mm ==3 or mm == 5 or mm == 7 or mm == 8 or mm == 10):
        return 0
    elif dd > 28 and mm == 2 and (yr%4)!=0:
        return 0
    elif (yr%4)==0 and mm == 2 and dd>29:
        return 0
    elif dd > 31 and mm == 12:
        return 0
    else:
        return 1
#*********************************************************************************************************************

