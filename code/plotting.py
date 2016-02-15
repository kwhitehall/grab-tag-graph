import glob
import os
import subprocess
import time
from datetime import timedelta, datetime

import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import numpy.ma as ma
from matplotlib import cm
from matplotlib.ticker import FuncFormatter, FormatStrFormatter
from netCDF4 import Dataset

import mccSearch


# **********************************************************************************************************************


def draw_graph(thisGraph, graphTitle, MAIN_DIRECTORY, edgeWeight=None):
    '''
    Purpose:: Utility function to draw graph in the hierachial format

    Input:: thisGraph: a Networkx directed graph representing a subgraph of connected cloud elements
            graphTitle: a string representing the graph title
            MAIN_DIRECTORY: a string representing the path to the main directory where the data generated is saved
            edgeWeight: (optional) a list of integers representing the edge weights in the graph

    Returns:: None

    Generates:: a plot of the thisGraph

    '''

    imgFilename = MAIN_DIRECTORY + '/images/' + graphTitle + '.gif'
    fig = plt.figure(facecolor='white', figsize=(16, 12))

    edge95 = [(u, v) for (u, v, d) in thisGraph.edges(data=True) if d['weight'] == edgeWeight[0]]
    edge90 = [(u, v) for (u, v, d) in thisGraph.edges(data=True) if d['weight'] == edgeWeight[1]]
    edegeOverlap = [(u, v) for (u, v, d) in thisGraph.edges(data=True) if d['weight'] == edgeWeight[2]]

    nx.write_dot(thisGraph, 'test.dot')
    plt.title(graphTitle)
    pos = nx.graphviz_layout(thisGraph, prog='dot')
    # draw graph in parts
    # nodes
    nx.draw_networkx_nodes(thisGraph, pos, with_labels=True, arrows=False)
    # edges
    nx.draw_networkx_edges(thisGraph, pos, edgelist=edge95, alpha=0.5, arrows=False)
    nx.draw_networkx_edges(thisGraph, pos, edgelist=edge90, edge_color='b', style='dashed', arrows=False)
    nx.draw_networkx_edges(thisGraph, pos, edgelist=edegeOverlap, edge_color='y', style='dashed', arrows=False)
    # labels
    nx.draw_networkx_labels(thisGraph, pos, arrows=False)
    plt.axis('off')
    plt.savefig(imgFilename, facecolor=fig.get_facecolor(), transparent=True)
    #do some clean up...and ensuring that we are in the right dir
    os.chdir((MAIN_DIRECTORY+'/'))
    subprocess.call('rm test.dot', shell=True)
# **********************************************************************************************************************


def display_size(finalMCCList, MAIN_DIRECTORY):
    '''
    Purpose:: To create a figure showing the area verse time for each MCS

    Inputs:: finalMCCList: a list of list of strings representing the list of nodes representing a MCC
            MAIN_DIRECTORY: a string representing the path to the main directory where the data generated is saved

    Returns:: None

    Generates:: a plot with the area contribution of each node in a feature in the list

    '''
    timeList = []
    count = 1
    imgFilename = ''
    minArea = 10000.0
    maxArea = 0.0
    eachNode = {}

    # for each node in the list, get the area information from the dictionary
    # in the graph and calculate the area

    if finalMCCList:
        for eachMCC in finalMCCList:
            # get the info from the node
            for node in eachMCC:
                eachNode = mccSearch.this_dict(node)
                timeList.append(eachNode['cloudElementTime'])

                if eachNode['cloudElementArea'] < minArea:
                    minArea = eachNode['cloudElementArea']
                if eachNode['cloudElementArea'] > maxArea:
                    maxArea = eachNode['cloudElementArea']

            #sort and remove duplicates
            timeList = list(set(timeList))
            timeList.sort()
            tdelta = timeList[1] - timeList[0]
            starttime = timeList[0] - tdelta
            endtime = timeList[-1] + tdelta
            timeList.insert(0, starttime)
            timeList.append(endtime)

            # plot info
            plt.close('all')
            title = 'Area distribution of the MCC over somewhere'
            fig = plt.figure(facecolor='white', figsize=(18, 10))
            fig, ax = plt.subplots(1, facecolor='white', figsize=(10, 10))

            # the data
            for node in eachMCC:
                eachNode = mccSearch.this_dict(node)
                if eachNode['cloudElementArea'] < 80000:
                    ax.plot(eachNode['cloudElementTime'], eachNode['cloudElementArea'], 'bo', markersize=10)
                elif eachNode['cloudElementArea'] >= 80000.00 and eachNode['cloudElementArea'] < 160000.00:
                    ax.plot(eachNode['cloudElementTime'], eachNode['cloudElementArea'], 'yo', markersize=20)
                else:
                    ax.plot(eachNode['cloudElementTime'], eachNode['cloudElementArea'], 'ro', markersize=30)

            #axes and labels
            maxArea += 1000.00
            ax.set_xlim(starttime, endtime)
            ax.set_ylim(minArea, maxArea)
            ax.set_ylabel('Area in km^2', fontsize=12)
            ax.set_title(title)
            ax.fmt_xdata = mdates.DateFormatter('%Y-%m-%d%H:%M:%S')
            fig.autofmt_xdate()

            plt.subplots_adjust(bottom=0.2)

            imgFilename = MAIN_DIRECTORY+'/images/' + str(count)+'MCS.gif'
            plt.savefig(imgFilename, facecolor=fig.get_facecolor(), transparent=True)

            # if time in not already in the time list, append it
            timeList = []
            count += 1
    return
# **********************************************************************************************************************


def display_precip(finalMCCList, MAIN_DIRECTORY):
    '''
    Purpose:: To create a figure showing the precip rate verse time for each MCS

    Inputs:: finalMCCList: a list of dictionaries representing a list of nodes representing a MCC
            MAIN_DIRECTORY: a string representing the path to the main directory where the data generated is saved

    Returns:: None

    Generates:: A plot with the precipitation distribution over each feature in the list

    '''
    timeList = []
    oriTimeList = []
    colorBarTime = []
    count = 1
    imgFilename = ''
    percentagePrecipitating = []  #0.0
    cloudElementArea = []
    nodes = []
    xy = []
    x = []
    y = []
    precip = []
    partialArea = []
    totalSize = 0.0
    firstTime = True
    xStart = 0.0
    yStart = 0.0

    # for each node in the list, get the area information from the dictionary
    # in the graph and calculate the area
    if finalMCCList:
        for eachMCC in finalMCCList:
            # Get the info from the node
            for node in eachMCC:
                eachNode = mccSearch.this_dict(node)
                if firstTime is True:
                    xStart = eachNode['cloudElementCenter'][1]  # lon
                    yStart = eachNode['cloudElementCenter'][0]  # lat
                timeList.append(eachNode['cloudElementTime'])
                percentagePrecipitating.append((eachNode['TRMMArea'] / eachNode['cloudElementArea']) * 100.0)
                cloudElementArea.append(eachNode['cloudElementArea'])
                nodes.append(eachNode['uniqueID'])
                x.append(eachNode['cloudElementCenter'][1])  # -xStart)
                y.append(eachNode['cloudElementCenter'][0])  # -yStart)

                firstTime = False

            #convert the timeList[] to list of floats
            for i in xrange(len(timeList)):
                colorBarTime.append(time.mktime(timeList[i].timetuple()))

            totalSize = sum(cloudElementArea)
            partialArea = [(a/totalSize)*30000 for a in cloudElementArea]

            # Plot info
            plt.close('all')

            title = 'Precipitation distribution of the MCS '
            fig, ax = plt.subplots(1, facecolor='white', figsize=(20, 7))

            cmap = cm.jet
            ax.scatter(x, y, s=partialArea, c=colorBarTime, edgecolors='none', marker='o', cmap=cmap)
            colorBarTime = []
            colorBarTime = list(set(timeList))
            colorBarTime.sort()
            cb = colorbar_index(ncolors=len(colorBarTime), nlabels=colorBarTime, cmap=cmap)

            # axes and labels
            ax.set_xlabel('Degrees Longtude', fontsize=12)
            ax.set_ylabel('Degrees Latitude', fontsize=12)
            ax.set_title(title)
            ax.grid(True)
            plt.subplots_adjust(bottom=0.2)

            for i, txt in enumerate(nodes):
                if cloudElementArea[i] >= 2400.00:
                    ax.annotate('%d' % percentagePrecipitating[i]+'%', (x[i], y[i]))
                precip = []

            imgFilename = MAIN_DIRECTORY+'/images/MCSprecip' + str(count)+'.gif'
            plt.savefig(imgFilename, facecolor=fig.get_facecolor(), transparent=True)

            # Reset for next image
            timeList = []
            percentagePrecipitating = []
            cloudElementArea = []
            x = []
            y = []
            colorBarTime = []
            nodes = []
            precip = []
            count += 1
            firstTime = True
    return
# **********************************************************************************************************************


def plot_accu_in_time_range(starttime, endtime, MAIN_DIRECTORY, tRes):
    '''
    Purpose:: Create accumulated precip plot within a time range given using all CEs

    Inputs:: starttime: a string representing the time to start the accumulations format yyyy-mm-dd_hh:mm:ss
             endtime: a string representing the time to end the accumulations format yyyy-mm-dd_hh:mm:ss
             MAIN_DIRECTORY: a string representing the path to the main directory where the data generated is saved
             tRes: a float representing the time res of the input data e.g. 30min=0.5
             
    Returns:: a netcdf file containing the accumulated precip for specified times a gif (generated in Grads)

    Generates:: a plot of the amount of precipitation accumulated between two dates for the domain

    TODO: pass of pick up from the NETCDF file the  lat, lon and resolution for generating the ctl file
    '''

    os.chdir((MAIN_DIRECTORY+'/TRMMnetcdfCEs/'))
    # Just incase the X11 server is giving problems
    subprocess.call('export DISPLAY=:0.0', shell=True)

    imgFilename = ''
    firstTime = True

    fileList = []
    sTime = datetime.strptime(starttime.replace('_', ' '), '%Y-%m-%d %H:%M:%S')
    eTime = datetime.strptime(endtime.replace('_', ' '), '%Y-%m-%d %H:%M:%S')
    thisTime = sTime

    while thisTime <= eTime:
        fileList = filter(os.path.isfile, glob.glob(('TRMM' + str(thisTime).replace(' ', '_') + '*' + '.nc')))
        for fname in fileList:
            cloudElementTRMMData = Dataset(fname, 'r', format='NETCDF4')
            precipRate = cloudElementTRMMData.variables['precipitation_Accumulation'][:]
            lats = cloudElementTRMMData.variables['latitude'][:]
            lons = cloudElementTRMMData.variables['longitude'][:]
            LONTRMM, LATTRMM = np.meshgrid(lons, lats)
            nygrdTRMM = len(LATTRMM[:, 0])
            nxgrdTRMM = len(LONTRMM[0, :])
            precipRate = ma.masked_array(precipRate, mask=(precipRate < 0.0))
            cloudElementTRMMData.close()

            if firstTime is True:
                accuPrecipRate = ma.zeros(precipRate.shape)
                firstTime = False

            accuPrecipRate += precipRate

        # Increment the time
        thisTime += timedelta(hours=tRes)

    # Create new netCDF file
    accuTRMMFile = MAIN_DIRECTORY+'/TRMMnetcdfCEs/accu'+starttime+'-'+endtime+'.nc'
    print 'accuTRMMFile ', accuTRMMFile
    # Write the file
    accuTRMMData = Dataset(accuTRMMFile, 'w', format='NETCDF4')
    accuTRMMData.description = 'Accumulated precipitation data'
    accuTRMMData.calendar = 'standard'
    accuTRMMData.conventions = 'COARDS'
    # Dimensions
    accuTRMMData.createDimension('time', None)
    accuTRMMData.createDimension('lat', nygrdTRMM)
    accuTRMMData.createDimension('lon', nxgrdTRMM)

    # Variables
    TRMMprecip = ('time', 'lat', 'lon',)
    times = accuTRMMData.createVariable('time', 'f8', ('time',))
    times.units = 'hours since ' + starttime[:-6]
    latitude = accuTRMMData.createVariable('latitude', 'f8', ('lat',))
    longitude = accuTRMMData.createVariable('longitude', 'f8', ('lon',))
    rainFallacc = accuTRMMData.createVariable('precipitation_Accumulation', 'f8', TRMMprecip)
    rainFallacc.units = 'mm'

    longitude[:] = LONTRMM[0, :]
    longitude.units = 'degrees_east'
    longitude.long_name = 'Longitude'

    latitude[:] = LATTRMM[:, 0]
    latitude.units = 'degrees_north'
    latitude.long_name = 'Latitude'

    rainFallacc[:] = accuPrecipRate[:]

    accuTRMMData.close()

    # Generate the image with GrADS
    # the ctl file
    subprocess.call('rm acc.ctl', shell=True)
    subprocess.call('touch acc.ctl', shell=True)
    replaceExpDset = 'echo DSET ' + accuTRMMFile + ' >> acc.ctl'
    subprocess.call(replaceExpDset, shell=True)
    subprocess.call('echo "OPTIONS yrev little_endian template" >> acc.ctl', shell=True)
    subprocess.call('echo "DTYPE netcdf" >> acc.ctl', shell=True)
    subprocess.call('echo "UNDEF  0" >> acc.ctl', shell=True)
    subprocess.call('echo "TITLE  TRMM MCS accumulated precipitation" >> acc.ctl', shell=True)
    replaceExpXDef = 'echo XDEF ' + str(nxgrdTRMM) + ' LINEAR ' + str(min(lons)) + ' ' + \
        str((max(lons)-min(lons))/nxgrdTRMM) + ' >> acc.ctl'
    subprocess.call(replaceExpXDef, shell=True)
    replaceExpYDef = 'echo YDEF '+str(nygrdTRMM)+' LINEAR '+str(min(lats)) + ' ' + \
        str((max(lats)-min(lats))/nygrdTRMM)+' >>acc.ctl'
    subprocess.call(replaceExpYDef, shell=True)
    #subprocess.call('echo "XDEF 384 LINEAR  -8.96875 0.036378335 " >> acc.ctl', shell=True)
    #subprocess.call('echo "YDEF 384 LINEAR 5.03515625 0.036378335 " >> acc.ctl', shell=True)
    subprocess.call('echo "ZDEF   01 LEVELS 1" >> acc.ctl', shell=True)
    subprocess.call('echo "TDEF 99999 linear 31aug2009 1hr" >> acc.ctl', shell=True)
    subprocess.call('echo "VARS 1" >> acc.ctl', shell=True)
    subprocess.call('echo "precipitation_Accumulation=>precipAcc     1  t,y,x    precipAccu" >> acc.ctl', shell=True)
    subprocess.call('echo "ENDVARS" >> acc.ctl', shell=True)
    # Generate GrADS script
    imgFilename = MAIN_DIRECTORY+'/images/accu'+starttime+'-'+endtime+'.gif'
    subprocess.call('rm accuTRMM1.gs', shell=True)
    subprocess.call('touch accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'reinit''\'" >> accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'open acc.ctl ''\'" >> accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'set grads off''\'" >> accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'set mpdset hires''\'" >> accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'set gxout shaded''\'" >> accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'set datawarn off''\'" >> accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'d precipacc''\'" >> accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'draw title TRMM Accumulated Precipitation [mm]''\'" >> accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'run cbarn''\'" >> accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'printim '+imgFilename + ' x1000 y800 white''\'" >> accuTRMM1.gs', shell=True)
    subprocess.call('echo "''\'quit''\'" >> accuTRMM1.gs', shell=True)
    gradscmd = 'grads -blc ' + '\'run accuTRMM1.gs''\''
    subprocess.call(gradscmd, shell=True)

    # Clean up
    subprocess.call('rm accuTRMM1.gs', shell=True)
    subprocess.call('rm acc.ctl', shell=True)

    return
# **********************************************************************************************************************


def plot_accu_TRMM(finalMCCList, MAIN_DIRECTORY):
    '''
    Purpose:: (1) generate a file with the accumulated precipitation for the MCS (2) generate the appropriate image
            
    Input:: finalMCCList: a list of dictionaries representing a list of nodes representing a MCC
            MAIN_DIRECTORY: a string representing the path to the main directory where the data generated is saved

    Returns:: a netcdf file containing the accumulated precip
    
    Generates: A plot (generated in Grads) of the accumulated precipitation

    '''

    os.chdir((MAIN_DIRECTORY+'/TRMMnetcdfCEs'))
    fname = ''
    imgFilename = ''
    firstPartName = ''
    firstTime = True
    replaceExpXDef = ''
    thisNode = ''

    subprocess.call('export DISPLAY=:0.0', shell=True)

    # Generate the file name using MCCTimes
    # if the file name exists, add it to the accTRMM file
    for path in finalMCCList:
        firstTime = True
        for eachNode in path:
            thisNode = mccSearch.this_dict(eachNode)
            fname = 'TRMM' + str(thisNode['cloudElementTime']).replace(' ', '_') + thisNode['uniqueID'] + '.nc'

            if os.path.isfile(fname):
                # Open NetCDF file add info to the accu
                cloudElementTRMMData = Dataset(fname, 'r', format='NETCDF4')
                precipRate = cloudElementTRMMData.variables['precipitation_Accumulation'][:]
                lats = cloudElementTRMMData.variables['latitude'][:]
                lons = cloudElementTRMMData.variables['longitude'][:]
                LONTRMM, LATTRMM = np.meshgrid(lons, lats)
                nygrdTRMM = len(LATTRMM[:, 0])
                nxgrdTRMM = len(LONTRMM[0, :])
                precipRate = ma.masked_array(precipRate, mask=(precipRate < 0.0))
                cloudElementTRMMData.close()

                if firstTime is True:
                    firstPartName = str(thisNode['uniqueID'])+str(thisNode['cloudElementTime']).replace(' ', '_')+'-'
                    accuPrecipRate = ma.zeros(precipRate.shape)
                    firstTime = False

                accuPrecipRate += precipRate

        imgFilename = MAIN_DIRECTORY+'/images/MCS_'+firstPartName+str(thisNode['cloudElementTime']).replace(' ', '_')+'.gif'
        # Create new netCDF file
        accuTRMMFile = MAIN_DIRECTORY+'/TRMMnetcdfCEs/accu'+firstPartName+str(thisNode['cloudElementTime']).replace(' ', '_')+'.nc'
        # Write the file
        accuTRMMData = Dataset(accuTRMMFile, 'w', format='NETCDF4')
        accuTRMMData.description = 'Accumulated precipitation data'
        accuTRMMData.calendar = 'standard'
        accuTRMMData.conventions = 'COARDS'
        # dimensions
        accuTRMMData.createDimension('time', None)
        accuTRMMData.createDimension('lat', nygrdTRMM)
        accuTRMMData.createDimension('lon', nxgrdTRMM)

        # variables
        TRMMprecip = ('time', 'lat', 'lon',)
        times = accuTRMMData.createVariable('time', 'f8', ('time',))
        times.units = 'hours since ' + str(thisNode['cloudElementTime']).replace(' ', '_')[:-6]
        latitude = accuTRMMData.createVariable('latitude', 'f8', ('lat',))
        longitude = accuTRMMData.createVariable('longitude', 'f8', ('lon',))
        rainFallacc = accuTRMMData.createVariable('precipitation_Accumulation', 'f8', TRMMprecip)
        rainFallacc.units = 'mm'

        longitude[:] = LONTRMM[0, :]
        longitude.units = 'degrees_east'
        longitude.long_name = 'Longitude'

        latitude[:] = LATTRMM[:, 0]
        latitude.units = 'degrees_north'
        latitude.long_name = 'Latitude'

        rainFallacc[:] = accuPrecipRate[:]

        accuTRMMData.close()

        # Generate the image with GrADS
        subprocess.call('rm acc.ctl', shell=True)
        subprocess.call('touch acc.ctl', shell=True)
        replaceExpDset = 'echo DSET ' + accuTRMMFile + ' >> acc.ctl'
        subprocess.call(replaceExpDset, shell=True)
        subprocess.call('echo "OPTIONS yrev little_endian template" >> acc.ctl', shell=True)
        subprocess.call('echo "DTYPE netcdf" >> acc.ctl', shell=True)
        subprocess.call('echo "UNDEF  0" >> acc.ctl', shell=True)
        subprocess.call('echo "TITLE  TRMM MCS accumulated precipitation" >> acc.ctl', shell=True)
        replaceExpXDef = 'echo XDEF ' + str(nxgrdTRMM) + ' LINEAR ' + str(min(lons)) + ' ' + \
            str((max(lons)-min(lons))/nxgrdTRMM) + ' >> acc.ctl'
        subprocess.call(replaceExpXDef, shell=True)
        replaceExpYDef = 'echo YDEF '+str(nygrdTRMM)+' LINEAR '+str(min(lats)) + ' ' + \
            str((max(lats)-min(lats))/nygrdTRMM)+' >>acc.ctl'
        subprocess.call(replaceExpYDef, shell=True)
        subprocess.call('echo "ZDEF   01 LEVELS 1" >> acc.ctl', shell=True)
        subprocess.call('echo "TDEF 99999 linear 31aug2009 1hr" >> acc.ctl', shell=True)
        subprocess.call('echo "VARS 1" >> acc.ctl', shell=True)
        subprocess.call('echo "precipitation_Accumulation=>precipAcc     1  t,y,x    precipAccu" >> acc.ctl', shell=True)
        subprocess.call('echo "ENDVARS" >> acc.ctl', shell=True)

        # Generate GrADS script
        subprocess.call('rm accuTRMM1.gs', shell=True)
        subprocess.call('touch accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'reinit''\'" >> accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'open acc.ctl ''\'" >> accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'set grads off''\'" >> accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'set mpdset hires''\'" >> accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'set gxout shaded''\'" >> accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'set datawarn off''\'" >> accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'d precipacc''\'" >> accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'draw title TRMM Accumulated Precipitation [mm]''\'" >> accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'run cbarn''\'" >> accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'printim '+imgFilename + ' x1000 y800 white''\'" >> accuTRMM1.gs', shell=True)
        subprocess.call('echo "''\'quit''\'" >> accuTRMM1.gs', shell=True)
        gradscmd = 'grads -blc ' + '\'run accuTRMM1.gs''\''
        subprocess.call(gradscmd, shell=True)

        # Clean up
        subprocess.call('rm accuTRMM1.gs', shell=True)
        subprocess.call('rm acc.ctl', shell=True)

    return
# **********************************************************************************************************************


def plot_histogram(precip, aTitle, imgFilename, xlabel, ylabel, numBins):
    '''
    Purpose:: To create plots (histograms) of the data entered in aList

    Input:: precip: a list of the values to generate the histogram
        aTitle: a string representing the title and the name of the plot e.g. 'Average area [km^2]'
        imgFilename: a string representing filename for the image
        xlabel: a string for the label on the x-axis
        ylabel: a string for the label on the y-axis
        numBins: an integer representing the num of bins to use in the histogram

    Returns:: None

    Generates:: A histogram plot of the precip array values passed

    '''

    fig, ax = plt.subplots(1, facecolor='white', figsize=(7, 5))

    n, binsdg = np.histogram(precip, numBins)
    wid = binsdg[1:] - binsdg[:-1]
    plt.bar(binsdg[:-1], n/float(len(precip)), width=wid)

    # Make percentage plot
    formatter = FuncFormatter(to_percent)
    plt.xlim(min(binsdg), max(binsdg))
    ax.set_xticks(binsdg)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(aTitle)
    # Set the formatter
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%0.0f'))
    plt.savefig(imgFilename, transparent=True)

    return
# **********************************************************************************************************************


def plot_precip_histograms(finalMCCList, MAIN_DIRECTORY):
    '''
    Purpose:: To create plots (histograms) of the each TRMMnetcdfCEs files

    Input:: finalMCCList: a list of dictionaries representing a list of nodes representing a MCC
        MAIN_DIRECTORY: a string representing the path to the main directory where the data generated is saved

    Returns:: None

    Generates:: Histogram plots for each feature in the list passed
        
    '''
    numBins = 5
    precip = []
    imgFilename = ' '
    startTime = ' '
    firstTime = True

    if finalMCCList:

        for eachMCC in finalMCCList:
            firstTime = True

            for node in eachMCC:
                eachNode = mccSearch.this_dict(node)
                thisTime = eachNode['cloudElementTime']

                thisFileName = MAIN_DIRECTORY+'/TRMMnetcdfCEs/TRMM' + str(thisTime).replace(' ', '_') + \
                    eachNode['uniqueID'] + '.nc'
                TRMMData = Dataset(thisFileName, 'r', format='NETCDF4')
                precipRate = TRMMData.variables['precipitation_Accumulation'][:, :, :]
                cloudElementPrecipRate = precipRate[0, :, :]
                TRMMData.close()

                if firstTime is True:
                    totalPrecip = np.zeros(cloudElementPrecipRate.shape)
                    startTime = str(thisTime) + eachNode['uniqueID']
                    firstTime = False

                totalPrecip = np.add(totalPrecip, precipRate)

                for _, value in np.ndenumerate(cloudElementPrecipRate):
                    if value != 0.0:
                        precip.append(value)

                plt.close('all')
                title = 'TRMM precipitation distribution for ' + str(thisTime)
                imgFilename = MAIN_DIRECTORY+'/images/'+str(thisTime)+eachNode['uniqueID']+'TRMMMCS.gif'
                xlabel = 'Precipitation [mm]'
                ylabel = 'Area [km^2]'
                numBins = 5
                if all(value == 0 for value in precip):
                    print 'NO Precipitation at '+str(thisTime)+eachNode['uniqueID']
                else:
                    plot_histogram(precip, title, imgFilename, xlabel, ylabel, numBins)

                precip = []

            # The image at the end of each MCS
            title = 'TRMM precipitation distribution for '+startTime+' to '+str(thisTime)
            imgFilename = MAIN_DIRECTORY+'/images/'+startTime+'_'+str(thisTime)+eachNode['uniqueID']+'TRMMMCS.gif'
            xlabel = 'Precipitation [mm]'
            ylabel = 'Area [km^2]'
            numBins = 10
            for _, value in np.ndenumerate(totalPrecip):
                if value != 0.0:
                    precip.append(value)

            if all(value == 0 for value in precip):
                print 'No precipitation for MCS starting at '+startTime+' and ending at '+str(thisTime)+eachNode['uniqueID']
            else:
                plot_histogram(precip, title, imgFilename, xlabel, ylabel, numBins)

            precip = []

    return
# **********************************************************************************************************************
#           PLOTTING UTIL SCRIPTS
# **********************************************************************************************************************


def to_percent(y, position):
    '''
    Purpose::
        Utility script for generating the y-axis for plots
    '''
    return str(100*y) + '%'
# **********************************************************************************************************************


def colorbar_index(ncolors, nlabels, cmap):
    '''
    Purpose::
        Utility script for creating a colorbar
        Taken from http://stackoverflow.com/questions/18704353/correcting-matplotlib-colorbar-ticks

    Inputs:: ncolors: an int representing the number of colors in the colorbar
            nlabels: an  int representing the number of labels to be assigned in the colorbar
            cmap: a matplotlib cmap object representing which cmap style to use
    
    Returns:: None

    '''

    cmap = cmap_discretize(cmap, ncolors)
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar(mappable)
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    colorbar.set_ticklabels(nlabels)
    return

# **********************************************************************************************************************


def cmap_discretize(cmap, N):
    '''
    Taken from: http://stackoverflow.com/questions/18704353/correcting-matplotlib-colorbar-ticks
    http://wiki.scipy.org/Cookbook/Matplotlib/ColormapTransformations
    Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    '''

    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i-1, ki], colors_rgba[i, ki]) for i in xrange(N+1)]
    # Return colormap object.
    return mcolors.LinearSegmentedColormap(cmap.name + '_%d' % N, cdict, 1024)
# **********************************************************************************************************************

