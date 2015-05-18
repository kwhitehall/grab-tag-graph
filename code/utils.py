import glob
import itertools
import os
import re
import subprocess

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from scipy.ndimage import map_coordinates


def do_regrid(inputGrid, lat, lon, lat2, lon2, order=1, mdi=-999999999):
    '''
     Perform regridding from one set of lat,lon values onto a new set (lat2,lon2)

     Input::
         inputGrid          - the variable to be regridded
         lat,lon    - original co-ordinates corresponding to inputGrid values
         lat2,lon2  - new set of latitudes and longitudes that you want to regrid inputGrid onto
         order      - (optional) interpolation order 1=bi-linear, 3=cubic spline
         mdi        - (optional) fill value for missing data (used in creation of masked array)

     Output::
         outputGrid  - inputGrid regridded onto the new set of lat2,lon2

    '''

    nlat = inputGrid.shape[0]
    nlon = inputGrid.shape[1]

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

    """
    Notes on dealing with MDI when regridding data.
      Method adopted here:

        *  Use bilinear interpolation of data by default (but user can specify
        other order using order=... in call)

        *  Perform bilinear interpolation of data, and of mask.

        *  To be conservative, new grid point which contained some missing data
        on the old grid is set to missing data.  This is achieved by looking
        for any non-zero interpolated mask values.

        *  To avoid issues with bilinear interpolation producing strong
        gradients leading into the MDI, set values at MDI points to mean data
        value so little gradient visible = not ideal, but acceptable for now.

    """
    # Set values in MDI so that similar to surroundings so don't produce large
    # gradients when interpolating
    # Preserve MDI mask, by only changing data part of masked array object.
    for shift in (-1, 1):
        for axis in (0, 1):
            qShifted = np.roll(inputGrid, shift=shift, axis=axis)
            idx = ~qShifted.mask * inputGrid.mask
            inputGrid.data[idx] = qShifted[idx]

    # Now we actually interpolate
    # map_coordinates does cubic interpolation by default,
    # use "order=1" to preform bilinear interpolation instead...
    outputGrid = map_coordinates(inputGrid, [lati, loni], order=order)
    outputGrid = outputGrid.reshape([nlat2, nlon2])

    # Set values to missing data outside of original domain
    outputGrid = ma.masked_array(outputGrid, mask=np.logical_or(np.logical_or(lat2 >= lat.max(),
                                                              lat2 <= lat.min()),
                                                np.logical_or(lon2 <= lon.min(),
                                                              lon2 >= lon.max())))

    # Make second map using nearest neighbour interpolation -use this to determine locations with MDI and mask these
    qmdi = np.zeros_like(inputGrid)
    qmdi[inputGrid.mask == True] = 1.
    qmdi[inputGrid.mask == False] = 0.
    qmdiR = map_coordinates(qmdi, [lati, loni], order=order)
    qmdiR = qmdiR.reshape([nlat2, nlon2])
    mdimask = (qmdiR != 0.0)

    # Combine missing data mask, with outside domain mask define above.
    outputGrid.mask = np.logical_or(mdimask, outputGrid.mask)

    return outputGrid
#******************************************************************
def find_nearest(thisArray, value):
    '''
    Purpose :: to determine the value within an array closes to
            another value

    Input ::
    Output::
    '''
    idx = (np.abs(thisArray-value)).argmin()
    return thisArray[idx]
#******************************************************************
def find_time(curryr, currmm, currdd, currhr):
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
            currmm += 1
            currdd = 1
        elif currdd > 31 and (currmm == 1 or currmm == 3 or currmm == 5 or currmm == 7 or currmm == 8 or currmm == 10):
            currmm += 1
            currdd = 1
        elif currdd > 31 and currmm == 12:
            currmm = 1
            currdd = 1
            curryr += 1
        elif currdd > 28 and currmm == 2 and (curryr%4) != 0:
            currmm = 3
            currdd = 1
        elif (curryr%4) == 0 and currmm == 2 and currdd > 29:
            currmm = 3
            currdd = 1

    if currmm < 10:
        currmmStr = "0"+str(currmm)
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
def preprocessing_merg(mergDirname):
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

    os.chdir((mergDirname+'/'))
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
        ftime = re.search('\_(.*)\_', fname).group(1)

        year = ftime[0:4]
        digitMonth = ftime[4:6]
        day = ftime[6:8]
        hour = ftime[8:10]

        digitMonthToStringMap = {'01' : 'Jan',
                                 '02' : 'Feb',
                                 '03' : 'Mar',
                                 '04' : 'Apr',
                                 '05' : 'May',
                                 '06' : 'Jun',
                                 '07' : 'Jul',
                                 '08' : 'Aug',
                                 '09' : 'Sep',
                                 '10' : 'Oct',
                                 '11' : 'Nov',
                                 '12' : 'Dec'
                                }

        month = digitMonthToStringMap[digitMonth]

        subprocess.call('rm merg.ctl', shell=True)
        subprocess.call('touch merg.ctl', shell=True)
        replaceExpDset = 'echo DSET ' + fname +' >> merg.ctl'
        replaceExpTdef = 'echo TDEF 99999 LINEAR '+hour+'z'+day+month+year +' 30mn' +' >> merg.ctl'
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
        #lats4D = 'lats4d -v -q -lat '+LATMIN + ' ' +LATMAX +' -lon ' +LONMIN +' ' +LONMAX +' -time '+hour+'Z'+day+month+year + ' -func @+75 ' + '-i merg.ctl' + ' -o ' + fname

        #lats4D = 'lats4d -v -q -lat -40 -15 -lon 10 40 -time '+hour+'Z'+day+month+year + ' -func @+75 ' + '-i merg.ctl' + ' -o ' + fname
        #lats4D = 'lats4d -v -q -lat -5 40 -lon -90 60 -func @+75 ' + '-i merg.ctl' + ' -o ' + fname

        #gradscmd = 'grads -blc ' + '\'' +lats4D + '\''
        #run grads and lats4d command
        #subprocess.call(gradscmd, shell=True)
        imgFilename = hour+'Z'+day+month+year+'.gif'
        temp_masked_images(imgFilename)

    #when all the files have benn converted, mv the netcdf files
    #subprocess.call('mkdir mergNETCDF', shell=True)
    #subprocess.call('mv *.nc mergNETCDF', shell=True)
    #mv all images
    subprocess.call('mkdir mergImgs', shell=True)
    subprocess.call('mv *.gif mergImgs', shell=True)
    return
#******************************************************************
def post_processing_netcdf(dataset, dirName=None):
    '''

    TODO: UPDATE TO PICK UP LIMITS FROM FILE FOR THE GRADS SCRIPTS
    TODO: Include a default if dirName isn't provided.  If dirName is left as None
          then the string concatinations for file pathing will throw a TypeError exception.
    TODO: Remove the MAINDIRECTORY variable since it is no longer defined.

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
    imgFilename = ''
    frameList = []
    fileList = []
    var = ''
    firstTime = True
    lineNum = 1
    #Just incase the X11 server is giving problems
    subprocess.call('export DISPLAY=:0.0', shell=True)

    prevFrameNum = 0

    if dataset == 1:
        var = 'ch4'
        ctlTitle = 'TITLE MCC search Output Grid: Time  lat lon'
        ctlLine = 'brightnesstemp=\>ch4     1  t,y,x    brightnesstemperature'
        origsFile = coreDir+"/cs1.gs"
        subprocessCall = 'touch '+origsFile
        subprocess.call(subprocessCall, shell=True)
        write_c3_grad_script(origsFile)
        gsFile = coreDir+"/cs2.gs"
        sologsFile = coreDir+"/mergeCE.gs"
        lineNum = 32
    elif dataset == 2:
        var = 'precipAcc'
        ctlTitle = 'TITLE  TRMM MCS accumulated precipitation search Output Grid: Time  lat lon '
        ctlLine = 'precipitation_Accumulation=\>precipAcc     1  t,y,x    precipAccu'
        origsFile = coreDir+"/cs3.gs"
        subprocessCall = 'touch '+origsFile
        subprocess.call(subprocessCall, shell=True)
        write_c1_grad_script(origsFile)
        gsFile = coreDir+"/cs4.gs"
        sologsFile = coreDir+"/TRMMCE.gs"
        lineNum = 10
    elif dataset == 3:
        var = 'ch4'
        ctlTitle = 'TITLE MERG DATA'
        ctlLine = 'ch4=\>ch4     1  t,y,x    brightnesstemperature'
        origsFile = coreDir+"/cs1.gs"
        subprocessCall = 'touch '+origsFile
        subprocess.call(subprocessCall, shell=True)
        write_c3_grad_script(origsFile)
        sologsFile = coreDir+"/infrared.gs"
        lineNum = 32


    #sort files
    os.chdir((dirName+'/'))
    try:
        os.makedirs('ctlFiles')
    except OSError:
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
            fileData = Dataset(fname, 'r', format='NETCDF4')
            lats = fileData.variables['latitude'][:]
            lons = fileData.variables['longitude'][:]
            lonData, latData = np.meshgrid(lons, lats)
            nygrd = len(latData[:, 0])
            nxgrd = len(lonData[0, :])
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
        lineToWrite = 'echo '+ctlLine+' >> '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite = 'echo ENDVARS >>  '+ctlFile1
        subprocess.call(lineToWrite, shell=True)
        lineToWrite = 'echo  >>  '+ctlFile1
        subprocess.call(lineToWrite, shell=True)

        #create plot of just that data
        subprocessCall = 'cp '+ origsFile+' '+sologsFile
        subprocess.call(subprocessCall, shell=True)

        imgFilename = fnameNoExtension + '.gif'

        displayCmd = '\''+'d '+ var+'\''+'\n'
        newFileCmd = '\''+'open '+ ctlFile1+'\''+'\n'
        colorbarCmd = '\''+'run cbarn'+'\''+'\n'
        printimCmd = '\''+'printim '+MAINDIRECTORY+'/images/'+imgFilename+' x800 y600 white\''+'\n'
        quitCmd = '\''+'quit'+'\''+'\n'

        gradsScript = open(sologsFile, 'r+')
        lines1 = gradsScript.readlines()
        gradsScript.seek(0)
        lines1.insert((1), newFileCmd)
        lines1.insert((lineNum+1), displayCmd)
        lines1.insert((lineNum+2), colorbarCmd)
        lines1.insert((lineNum+3), printimCmd)
        lines1.insert((lineNum + 4), quitCmd)
        gradsScript.writelines(lines1)
        gradsScript.close()
        #run the script
        runGrads = 'run '+ sologsFile
        gradscmd = 'grads -blc ' + '\'' +runGrads + '\''+'\n'
        subprocess.call(gradscmd, shell=True)

        if dataset == 1 or dataset == 2:

            #TODO: for either dataset 1 or 2, write write_c3_grad_script and then use the for loop to gen the line to add at the end to display
            #at the end, run the GrADS script

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
                        ceNumber = fnameNoExtension.split('CE')[1]
                        imgFilename = fnameNoExtension.split('CE')[0] + '.gif'
                        ctlFile1 = dirName+'/ctlFiles/'+fnameNoExtension + '.ctl'

                        #build cs.gs script will all the CE ctl files and the appropriate display command
                        newVar = var+'.'+ceNumber
                        newDisplayCmd = '\''+'d '+ newVar+'\''+'\n'
                        newFileCmd = '\''+'open '+ ctlFile1+'\''+'\n'
                        gradsScript = open(gsFile, 'r+')
                        lines1 = gradsScript.readlines()
                        gradsScript.seek(0)
                        lines1.insert((1+count), newFileCmd)
                        lines1.insert((lineNum+count+1), newDisplayCmd)
                        gradsScript.writelines(lines1)
                        gradsScript.close()
                    count += 1

                colorbarCmd = '\''+'run cbarn'+'\''+'\n'
                printimCmd = '\''+'printim '+MAINDIRECTORY+'/images/'+imgFilename+' x800 y600 white\''+'\n'
                quitCmd = '\''+'quit'+'\''+'\n'
                gradsScript = open(gsFile, 'r+')
                lines1 = gradsScript.readlines()
                gradsScript.seek(0)
                lines1.insert((lineNum + (count*2)+1), colorbarCmd)
                lines1.insert((lineNum + (count*2)+2), printimCmd)
                lines1.insert((lineNum + (count*2)+3), quitCmd)
                gradsScript.writelines(lines1)
                gradsScript.close()

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
def temp_masked_images(imgFilename):
    '''
    Purpose::
        To generate temperature-masked images for a first pass verification

    Input::
        imgFilename: filename for the gif file

    Output::
        None - Gif images for each file of T_bb less than 250K are generated in folder called mergImgs

    Assumptions::
       Same as for preprocessing_merg
       1 GrADS (http://www.iges.org/grads/gadoc/) and lats4D (http://opengrads.org/doc/scripts/lats4d/)
         have been installed on the system and the user can access
       2 User can write files in location where script is being called
       3 the files havent been unzipped
    '''

    subprocess.call('rm temp_masked_images.gs', shell=True)
    subprocess.call('touch temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'open merg.ctl''\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'set mpdset hires''\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'set lat -5 30''\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'set lon -40 30''\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'set cint 10''\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'set clevs 190 200 210 220 230 240 250''\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'set gxout shaded''\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'d ch4+75''\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'run cbarn''\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'draw title Masked Temp @ '+imgFilename +'\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'printim '+imgFilename +' x1000 y800''\'" >> temp_masked_images.gs', shell=True)
    subprocess.call('echo "''\'quit''\'" >> temp_masked_images.gs', shell=True)
    gradscmd = 'grads -blc ' + '\'run temp_masked_images.gs''\''
    subprocess.call(gradscmd, shell=True)
    return
#******************************************************************
def valid_date(dataString):
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
    elif (dd < 0 or dd > 30) and (mm == 4 or mm == 6 or mm == 9 or mm == 11):
        return 0
    elif (dd < 0 or dd > 31) and (mm == 1 or mm == 3 or mm == 5 or mm == 7 or mm == 8 or mm == 10):
        return 0
    elif dd > 28 and mm == 2 and (yr%4) != 0:
        return 0
    elif (yr%4) == 0 and mm == 2 and dd > 29:
        return 0
    elif dd > 31 and mm == 12:
        return 0
    else:
        return 1
#*********************************************************************************************************************
def write_c3_grad_script(origsFile):
    '''
    Input:: a string representing the filename with full path to the GrADS script being created

    Output::

    Assumptions::

    '''
    subprocess.call('echo "''\'reinit''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set grads off''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set mpdset hires''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set gxout shaded''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set csmooth on''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set gxout shaded on''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set datawarn off''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set gxout shaded on''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set csmooth on''\'" >> '+origsFile, shell=True)
    return
#*********************************************************************************************************************
def write_c1_grad_script(origsFile):
    '''
    Input:: a string representing the filename with full path to the GrADS script being created

    Output::

    Assumptions::

    '''
    subprocess.call('echo "''\'reinit''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set grads off''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set mpdset hires''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set gxout shaded''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set csmooth on''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 16 255 255 255''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 17 108 20 156''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 18 77 50 183''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 19 48 83 213''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 20 22 107 236''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 21 0 193 254''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 22 42 166 255''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 23 66 197 249''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 24 92 226 255''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 25 124 255 249''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 26 132 252 204''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 27 135 252 145''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 28 151 255 130''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 29 209 255 128''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 30 255 246 117''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 31 255 189 58''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 32 249 136 6''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 33 241 110 0''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 34 212 93 1''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 35 208 68 0''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 36 182 48 10''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 37 163 29 2''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 38 138 15 0''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set rgb 39 255 255 255''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'set_plot(198,312,5)''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\' ''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\' ''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\' ''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'function set_plot(min,max,int)''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'    value = min''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'    cval=16''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'    c_levs = ''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'    c_cols = ''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'    while( value <= max )''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'      c_levs = c_levs ' ' value''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'      c_cols = c_cols ' ' cval''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'      value = value + int''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'      cval=cval+1''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'    endwhile''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'    c_cols=c_cols' 'cval-1''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'    set c_levs = \'c_levs''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'    set c_cols = \'c_cols''\'" >> '+origsFile, shell=True)
    subprocess.call('echo "''\'return''\'" >> '+origsFile, shell=True)
#*********************************************************************************************************************


