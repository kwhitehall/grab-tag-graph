import datetime
from datetime import timedelta, datetime

import numpy as np

import mccSearch
#**********************************************************************************************************************
#
#             METRICS FUNCTIONS FOR MERG.PY
# TODO: rewrite these metrics so that they use the data from the
#   file instead of the graph as this reduce mem resources needed
#
#
#**********************************************************************************************************************
def average_duration(allMCCtimes):
    '''
    Purpose:: To determine the average MCC length for the period

    Inputs:: allMCCtimes {MCCtimes, starttime, endtime, duration): a list of dictionaries
            of MCC temporal details for each MCC in the period considered

    Returns:: a floating-point representing the average duration of a MCC in the period

    Assumptions::

    '''

    return sum([MCC['duration'] for MCC in allMCCtimes], timedelta(seconds=0))/len(allMCCtimes)
#**********************************************************************************************************************
def average_feature_size(finalMCCList):
    '''
    Purpose:: To determine the average MCC size for the period

    Inputs:: finalMCCList: a list of list of strings representing  a list of list of nodes representing a MCC

    Returns::a floating-point representing the average area of a MCC in the period

    Assumptions::

    '''
    thisMCC = 0.0
    thisMCCAvg = 0.0

    #for each node in the list, get the are information from the dictionary
    #in the graph and calculate the area
    for eachPath in finalMCCList:
        for eachNode in eachPath:
            thisMCC += mccSearch.this_dict(eachNode)['cloudElementArea']

        thisMCCAvg += (thisMCC/len(eachPath))
        thisMCC = 0.0

    #calcuate final average
    return thisMCCAvg/(len(finalMCCList))
#**********************************************************************************************************************
def average_time(allTimes):
    '''
    Purpose:: To determine the average time in a list of datetimes e.g. of use is finding avg starttime,
    
    Inputs:: allTimes: a list of datetimes representing all of a given event e.g. start time

    Returns:: a floating-point number representing the average of the times given

    '''
    avgTime = 0

    for aTime in allTimes:
        avgTime += aTime.second + 60*aTime.minute + 3600*aTime.hour

    if len(allTimes) > 1:
        avgTime /= len(allTimes)

    rez = str(avgTime/3600) + ' ' + str((avgTime%3600)/60) + ' ' + str(avgTime%60)
    return datetime.strptime(rez, '%H %M %S')
#**********************************************************************************************************************
def common_feature_size(finalMCCList):
    '''
    Purpose:: To determine the common (mode) MCC size for the period

    Inputs:: finalMCCList: a list of list of strings representing the list of nodes representing a MCC

    Returns:: hist: an array representing the values of the histogram
              bin_edges: an array of floating-point representing the bin edges for the histogram

    Assumptions::

    '''
    thisMCC = 0.0
    thisMCCAvg = []

    #for each node in the list, get the area information from the dictionary
    #in the graph and calculate the area
    for eachPath in finalMCCList:
        for eachNode in eachPath:
            thisMCC += eachNode['cloudElementArea']

        thisMCCAvg.append(thisMCC/len(eachPath))
        thisMCC = 0.0

    #calcuate
    hist, bin_edges = np.histogram(thisMCCAvg)
    return hist, bin_edges
#**********************************************************************************************************************
def create_text_file(finalMCCList, identifier, MAIN_DIRECTORY, OUTER_CLOUD_SHIELD_AREA, TRES):
    '''
    Purpose::
        Create a text file with information about the MCS
        This function is expected to be especially of use regarding long term record checks

    Inputs::
        finalMCCList: a list of dictionaries representing a list of nodes representing a MCC
        identifier: an integer representing the type of list that has been entered...this is for creating file purposes
            1 - MCCList; 2- mcsList

    Returns:: None

    Generates::
        mcsUserFile: a user readable text file with all information about each MCS
        mcsSummaryFile: a user readable text file with the summary of the MCS
        mcsPostFile: a user readable text file with each MCS identified

    Assumptions::

    TODOs: provide visualizations for the calculations used to create the textFile
    '''

    durations = 0.0
    startTimes = []
    endTimes = []
    averagePropagationSpeed = 0.0
    speedCounter = 0
    maxArea = 0.0
    amax = 0.0
    avgMaxArea = []
    maxAreaCounter = 0.0
    maxAreaTime = ''
    eccentricity = 0.0
    firstTime = True
    matureFlag = True
    timeMCSMatures = ''
    maxCEprecipRate = 0.0
    minCEprecipRate = 0.0
    averageArea = 0.0
    averageAreaCounter = 0
    durationOfMatureMCC = 0
    avgMaxPrecipRate = 0.0
    avgMaxPrecipRateCounter = 0
    avgMinPrecipRate = 0.0
    avgMinPrecipRateCounter = 0
    cloudElementSpeed = 0.0
    mcsSpeed = 0.0
    mcsSpeedCounter = 0
    mcsPrecipTotal = 0.0
    avgmcsPrecipTotalCounter = 0
    bigPtotal = 0.0
    bigPtotalCounter = 0
    allPropagationSpeeds = []
    averageAreas = []
    areaAvg = 0.0
    avgPrecipTotal = 0.0
    avgPrecipTotalCounter = 0
    avgMaxmcsPrecipRate = 0.0
    avgMaxmcsPrecipRateCounter = 0
    avgMinmcsPrecipRate = 0.0
    avgMinmcsPrecipRateCounter = 0
    avgPrecipArea = []
    location = []
    avgPrecipAreaPercent = 0.0
    precipArea = 0.0
    precipAreaPercent = 0.0
    precipPercent = []
    precipCounter = 0
    precipAreaAvg = 0.0
    minSpeed = 0.0
    maxSpeed = 0.0

    if identifier == 1:
        mcsUserFile = open((MAIN_DIRECTORY + '/textFiles/MCCsUserFile.txt'), 'wb')
        mcsSummaryFile = open((MAIN_DIRECTORY + '/textFiles/MCCSummary.txt'), 'wb')
        mcsPostFile = open((MAIN_DIRECTORY + '/textFiles/MCCPostPrecessing.txt'), 'wb')

    if identifier == 2:
        mcsUserFile = open((MAIN_DIRECTORY + '/textFiles/MCSsUserFile.txt'), 'wb')
        mcsSummaryFile = open((MAIN_DIRECTORY + '/textFiles/MCSSummary.txt'), 'wb')
        mcsPostFile = open((MAIN_DIRECTORY + '/textFiles/MCSPostPrecessing.txt'), 'wb')

    for eachPath in finalMCCList:
        eachPath.sort(key=lambda nodeID: (len(nodeID.split('C')[0]), nodeID.split('C')[0], nodeID.split('CE')[1]))
        mcsPostFile.write('\n %s' % eachPath)

        startTime = mccSearch.this_dict(eachPath[0])['cloudElementTime']
        endTime = mccSearch.this_dict(eachPath[-1])['cloudElementTime']
        duration = (endTime - startTime) + timedelta(hours=TRES)

        # convert datatime duration to seconds and add to the total for the average duration of all MCS in finalMCCList
        durations += (duration.total_seconds())

        #durations += duration
        startTimes.append(startTime)
        endTimes.append(endTime)

        #get the precip info
        for eachNode in eachPath:
            thisNode = mccSearch.this_dict(eachNode)
            #set first time min 'fake' values
            if firstTime == True:
                minCEprecipRate = thisNode['CETRMMmin']
                avgMinmcsPrecipRate += thisNode['CETRMMmin']
                firstTime = False
            #calculate the speed
            # if thisNode['cloudElementArea'] >= OUTER_CLOUD_SHIELD_AREA:
            #   averagePropagationSpeed += find_cloud_element_speed(eachNode, eachPath)
            #   speedCounter +=1

            #Amax: find max area
            if thisNode['cloudElementArea'] > maxArea:
                maxArea = thisNode['cloudElementArea']
                maxAreaTime = str(thisNode['cloudElementTime'])
                eccentricity = thisNode['cloudElementEccentricity']
                location = thisNode['cloudElementCenter']

                #determine the time the feature matures
                if matureFlag == True:
                    timeMCSMatures = str(thisNode['cloudElementTime'])
                    matureFlag = False

            #find min and max precip rate
            if thisNode['CETRMMmin'] < minCEprecipRate:
                minCEprecipRate = thisNode['CETRMMmin']

            if thisNode['CETRMMmax'] > maxCEprecipRate:
                maxCEprecipRate = thisNode['CETRMMmax']

            #if MCC, then calculate for only mature phase else, calculate for full MCS
            if identifier == 1:
                #calculations for only the mature stage
                try:
                    if thisNode['nodeMCSIdentifier'] == 'M':
                        #calculate average area of the maturity feature only
                        averageArea += thisNode['cloudElementArea']
                        averageAreaCounter += 1
                        durationOfMatureMCC += 1
                        avgMaxPrecipRate += thisNode['CETRMMmax']
                        avgMaxPrecipRateCounter += 1
                        avgMinPrecipRate += thisNode['CETRMMmin']
                        avgMinPrecipRateCounter += 1
                        avgMaxmcsPrecipRate += thisNode['CETRMMmax']
                        avgMaxmcsPrecipRateCounter += 1
                        avgMinmcsPrecipRate += thisNode['CETRMMmin']
                        avgMinmcsPrecipRateCounter += 1

                        #the precip percentage (TRMM area/CE area)
                        if thisNode['cloudElementArea'] >= 0.0 and thisNode['TRMMArea'] >= 0.0:
                            precipArea += thisNode['TRMMArea']
                            avgPrecipArea.append(thisNode['TRMMArea'])
                            avgPrecipAreaPercent += (thisNode['TRMMArea']/thisNode['cloudElementArea'])
                            precipPercent.append((thisNode['TRMMArea']/thisNode['cloudElementArea']))
                            precipCounter += 1

                        #system speed for only mature stage
                        # cloudElementSpeed = find_cloud_element_speed(eachNode,eachPath)
                        # if cloudElementSpeed > 0.0 :
                        #   mcsSpeed += cloudElementSpeed
                        #   mcsSpeedCounter += 1
                except:
                    print 'MCS node has no nodeMCSIdentifier ', thisNode['uniqueID']
                    avgMaxmcsPrecipRate += thisNode['CETRMMmax']
                    avgMaxmcsPrecipRateCounter += 1
                    avgMinmcsPrecipRate += thisNode['CETRMMmin']
                    avgMinmcsPrecipRateCounter += 1
            else:
                averageArea += thisNode['cloudElementArea']
                averageAreaCounter += 1
                durationOfMatureMCC += 1
                avgMaxPrecipRate += thisNode['CETRMMmax']
                avgMaxPrecipRateCounter += 1
                avgMinPrecipRate += thisNode['CETRMMmin']
                avgMinPrecipRateCounter += 1
                avgMaxmcsPrecipRate += thisNode['CETRMMmax']
                avgMaxmcsPrecipRateCounter += 1
                avgMinmcsPrecipRate += thisNode['CETRMMmin']
                avgMinmcsPrecipRateCounter += 1

            #the precip percentage (TRMM area/CE area)
            if thisNode['cloudElementArea'] >= 0.0 and thisNode['TRMMArea'] >= 0.0:
                precipArea += thisNode['TRMMArea']
                avgPrecipArea.append(thisNode['TRMMArea'])
                avgPrecipAreaPercent += (thisNode['TRMMArea']/thisNode['cloudElementArea'])
                precipPercent.append((thisNode['TRMMArea']/thisNode['cloudElementArea']))
                precipCounter += 1

            #system speed for only mature stage
            # cloudElementSpeed = find_cloud_element_speed(eachNode,eachPath)
            # if cloudElementSpeed > 0.0 :
            #   mcsSpeed += cloudElementSpeed
            #   mcsSpeedCounter += 1

            #find accumulated precip
            if thisNode['cloudElementPrecipTotal'] > 0.0:
                mcsPrecipTotal += thisNode['cloudElementPrecipTotal']
                avgmcsPrecipTotalCounter += 1

        #A: calculate the average Area of the (mature) MCS
        if averageAreaCounter > 0: # and averageAreaCounter > 0:
            averageArea /= averageAreaCounter
            averageAreas.append(averageArea)

        #v: MCS speed
        if mcsSpeedCounter > 0: # and mcsSpeed > 0.0:
            mcsSpeed /= mcsSpeedCounter

        #smallP_max: calculate the average max precip rate (mm/h)
        if avgMaxmcsPrecipRateCounter > 0: #and avgMaxPrecipRate > 0.0:
            avgMaxmcsPrecipRate /= avgMaxmcsPrecipRateCounter

        #smallP_min: calculate the average min precip rate (mm/h)
        if avgMinmcsPrecipRateCounter > 0: #and avgMinPrecipRate > 0.0:
            avgMinmcsPrecipRate /= avgMinmcsPrecipRateCounter

        #smallP_avg: calculate the average precipitation (mm hr-1)
        if mcsPrecipTotal > 0.0: # and avgmcsPrecipTotalCounter> 0:
            avgmcsPrecipTotal = mcsPrecipTotal/avgmcsPrecipTotalCounter
            avgPrecipTotal += avgmcsPrecipTotal
            avgPrecipTotalCounter += 1

        #smallP_total = mcsPrecipTotal
        #precip over the MCS lifetime prep for bigP_total
        if mcsPrecipTotal > 0.0:
            bigPtotal += mcsPrecipTotal
            bigPtotalCounter += 1

        if maxArea > 0.0:
            avgMaxArea.append(maxArea)
            maxAreaCounter += 1

        #verage precipate area precentage (TRMM/CE area)
        if precipCounter > 0:
            avgPrecipAreaPercent /= precipCounter
            precipArea /= precipCounter


        #write stuff to file
        mcsUserFile.write('\n\n\nStarttime is: %s ' %(str(startTime)))
        mcsUserFile.write('\nEndtime is: %s ' %(str(endTime)))
        mcsUserFile.write('\nLife duration is %s hrs' %(str(duration)))
        mcsUserFile.write('\nTime of maturity is %s ' %(timeMCSMatures))
        mcsUserFile.write('\nDuration mature stage is: %s ' %durationOfMatureMCC*TRES)
        mcsUserFile.write('\nAverage area is: %.4f km^2 ' %(averageArea))
        mcsUserFile.write('\nMax area is: %.4f km^2 ' %(maxArea))
        mcsUserFile.write('\nMax area time is: %s ' %(maxAreaTime))
        mcsUserFile.write('\nEccentricity at max area is: %.4f ' %(eccentricity))
        mcsUserFile.write('\nCenter (lat,lon) at max area is: %.2f\t%.2f' %(location[0], location[1]))
        mcsUserFile.write('\nPropagation speed is %.4f ' %(mcsSpeed))
        mcsUserFile.write('\nMCS minimum preicip rate is %.4f mmh^-1' %(minCEprecipRate))
        mcsUserFile.write('\nMCS maximum preicip rate is %.4f mmh^-1' %(maxCEprecipRate))
        mcsUserFile.write('\nTotal precipitation during MCS is %.4f mm/lifetime' %(mcsPrecipTotal))
        mcsUserFile.write('\nAverage MCS precipitation is %.4f mm' %(avgmcsPrecipTotal))
        mcsUserFile.write('\nAverage MCS maximum precipitation is %.4f mmh^-1' %(avgMaxmcsPrecipRate))
        mcsUserFile.write('\nAverage MCS minimum precipitation is %.4f mmh^-1' %(avgMinmcsPrecipRate))
        mcsUserFile.write('\nAverage precipitation area is %.4f km^2 ' %(precipArea))
        mcsUserFile.write('\nPrecipitation area percentage of mature system %.4f percent ' %(avgPrecipAreaPercent*100))


        #append stuff to lists for the summary file
        if mcsSpeed > 0.0:
            allPropagationSpeeds.append(mcsSpeed)
            averagePropagationSpeed += mcsSpeed
            speedCounter += 1

        #reset vars for next MCS in list
        averageArea = 0.0
        averageAreaCounter = 0
        durationOfMatureMCC = 0
        mcsSpeed = 0.0
        mcsSpeedCounter = 0
        mcsPrecipTotal = 0.0
        avgMaxmcsPrecipRate = 0.0
        avgMaxmcsPrecipRateCounter = 0
        avgMinmcsPrecipRate = 0.0
        avgMinmcsPrecipRateCounter = 0
        firstTime = True
        matureFlag = True
        avgmcsPrecipTotalCounter = 0
        avgPrecipAreaPercent = 0.0
        precipArea = 0.0
        precipCounter = 0
        maxArea = 0.0
        maxAreaTime = ''
        eccentricity = 0.0
        timeMCSMatures = ''
        maxCEprecipRate = 0.0
        minCEprecipRate = 0.0
        location = []

    #LD: average duration
    if len(finalMCCList) > 1:
        durations /= len(finalMCCList)
        durations /= 3600.0 #convert to hours

        #A: average area
        areaAvg = sum(averageAreas)/ len(finalMCCList)
    #create histogram plot here
    #if len(averageAreas) > 1:
    #   plotHistogram(averageAreas, 'Average Area [km^2]', 'Area [km^2]')

    #Amax: average maximum area
    if maxAreaCounter > 0.0: #and avgMaxArea > 0.0 :
        amax = sum(avgMaxArea)/ maxAreaCounter
        #create histogram plot here
        #if len(avgMaxArea) > 1:
        #   plotHistogram(avgMaxArea, 'Maximum Area [km^2]', 'Area [km^2]')

    #v_avg: calculate the average propagation speed
    if speedCounter > 0:  # and averagePropagationSpeed > 0.0
        averagePropagationSpeed /= speedCounter

    #bigP_min: calculate the min rate in mature system
    if avgMinPrecipRate > 0.0: # and avgMinPrecipRateCounter > 0.0:
        avgMinPrecipRate /= avgMinPrecipRateCounter

    #bigP_max: calculate the max rate in mature system
    if avgMinPrecipRateCounter > 0.0: #and avgMaxPrecipRate >  0.0:
        avgMaxPrecipRate /= avgMaxPrecipRateCounter

    #bigP_avg: average total preicip rate mm/hr
    if avgPrecipTotalCounter > 0.0: # and avgPrecipTotal > 0.0:
        avgPrecipTotal /= avgPrecipTotalCounter

    #bigP_total: total precip rate mm/LD
    if bigPtotalCounter > 0.0: #and bigPtotal > 0.0:
        bigPtotal /= bigPtotalCounter

    #precipitation area percentage
    if len(precipPercent) > 0:
        precipAreaPercent = (sum(precipPercent)/len(precipPercent))*100.0

    #average precipitation area
    if len(avgPrecipArea) > 0:
        precipAreaAvg = sum(avgPrecipArea)/len(avgPrecipArea)
        #if len(avgPrecipArea) > 1:
        #   plotHistogram(avgPrecipArea, 'Average Rainfall Area [km^2]', 'Area [km^2]')


    sTime = str(average_time(startTimes))
    eTime = str(average_time(endTimes))
    if len(allPropagationSpeeds) > 1:
        maxSpeed = max(allPropagationSpeeds)
        minSpeed = min(allPropagationSpeeds)

    #write stuff to the summary file
    mcsSummaryFile.write('\nNumber of features is %d ' %(len(finalMCCList)))
    mcsSummaryFile.write('\nAverage duration is %.4f hrs ' %(durations))
    mcsSummaryFile.write('\nAverage startTime is %s ' %(sTime[-8:]))
    mcsSummaryFile.write('\nAverage endTime is %s ' %(eTime[-8:]))
    mcsSummaryFile.write('\nAverage size is %.4f km^2 ' %(areaAvg))
    mcsSummaryFile.write('\nAverage precipitation area is %.4f km^2 ' %(precipAreaAvg))
    mcsSummaryFile.write('\nAverage maximum size is %.4f km^2 ' %(amax))
    mcsSummaryFile.write('\nAverage propagation speed is %.4f ms^-1' %(averagePropagationSpeed))
    mcsSummaryFile.write('\nMaximum propagation speed is %.4f ms^-1 ' %(maxSpeed))
    mcsSummaryFile.write('\nMinimum propagation speed is %.4f ms^-1 ' %(minSpeed))
    mcsSummaryFile.write('\nAverage minimum precipitation rate is %.4f mmh^-1' %(avgMinPrecipRate))
    mcsSummaryFile.write('\nAverage maximum precipitation rate is %.4f mm h^-1' %(avgMaxPrecipRate))
    mcsSummaryFile.write('\nAverage precipitation is %.4f mm ' %(avgPrecipTotal))
    mcsSummaryFile.write('\nAverage total precipitation during MCSs is %.4f mm/LD ' %(bigPtotal))
    mcsSummaryFile.write('\nAverage precipitation area percentage is %.4f percent ' %(precipAreaPercent))

    mcsUserFile.close
    mcsSummaryFile.close
    mcsPostFile.close
    return
#**********************************************************************************************************************
def find_cloud_element_speed(node, mcsList, theList):
    '''
    Purpose:: To determine the speed of the cloud elements uses vector displacement deltaLat/deltaLon (y/x)

    Inputs:: node: a string representing the cloud element
        mcsList: a list of strings representing the feature

    Returns:: cloudElementSpeed: a floating-point number representing the speed of the cloud element

    '''

    deltaLon = 0.0
    deltaLat = 0.0
    cloudElementSpeed = []
    theSpeed = 0.0

    nodeLatLon = mccSearch.this_dict(node)['cloudElementCenter']

    for aNode in theList:
        if aNode in mcsList:
            #if aNode is part of the mcsList then determine distance
            aNodeLatLon = mccSearch.this_dict(aNode)['cloudElementCenter']
            #calculate CE speed
            #checking the lats
            # nodeLatLon[0] += 90.0
            # aNodeLatLon[0] += 90.0
            # deltaLat = (nodeLatLon[0] - aNodeLatLon[0])
            deltaLat = ((mccSearch.this_dict(node)['cloudElementCenter'][0] +90.0) - \
                (mccSearch.this_dict(aNode)['cloudElementCenter'][0]+90.0))
            # nodeLatLon[1] += 360.0
            # aNodeLatLon[1] += 360.0
            # deltaLon = (nodeLatLon[1] - aNodeLatLon[1])
            deltaLon = ((mccSearch.this_dict(node)['cloudElementCenter'][1]+360.0) - \
                (mccSearch.this_dict(aNode)['cloudElementCenter'][1]+360.0))
            #failsafe for movement only in one dir
            if deltaLat == 0.0:
                deltaLat = 1.0

            if deltaLon == 0.0:
                deltaLon = 1.0

            try:
                theSpeed = abs((((deltaLat/deltaLon)*LAT_DISTANCE*1000)/(TRES*3600))) #convert to s --> m/s
            except:
                theSpeed = 0.0

            cloudElementSpeed.append(theSpeed)

    if not cloudElementSpeed:
        return 0.0
    else:
        return min(cloudElementSpeed)
#**********************************************************************************************************************
def longest_duration(allMCCtimes):
    '''
    Purpose:: To determine the longest MCC for the period

    Inputs:: allMCCtimes: a list of dictionaries {MCCtimes, starttime, endtime, duration, area} representing
        all the detected MCC temporal details for each MCC in the period considered

    Returns:: a floating-point number representing the maximum duration MCC detected

    Assumptions::

    '''

    return max([MCC['duration'] for MCC in allMCCtimes])
#**********************************************************************************************************************
def number_of_features(finalMCCList):
    '''
    Purpose:: To count the number of MCCs found for the period

    Inputs:: finalMCCList: a list of list of strings representing a list of list of nodes in each MCC

    Returns:: an integer representing the number of MCCs found

    '''
    return len(finalMCCList)
#**********************************************************************************************************************
def precip_max_min(finalMCCList):
    '''
    TODO: this doesnt work the np.min/max function seems to be not working with the nonzero option
        ..possibly a problem upstream with cloudElementLatLonTRMM

    Purpose:: Precipitation maximum and min rates associated with each cloud elements in MCS
    
    Inputs:: finalMCCList: a list of dictionaries representing a list of nodes representing a MCC

    Returns:: mcsPrecip: a list representing the max and min rate for each cloud elements identified

    '''
    maxCEprecip = 0.0
    minCEprecip = 0.0
    mcsPrecip = []
    allmcsPrecip = []

    if finalMCCList:
        if type(finalMCCList[0]) is str: # len(finalMCCList) == 1:
            for node in finalMCCList:
                eachNode = mccSearch.this_dict(node)
                #cloudElementTRMM = eachNode['cloudElementLatLonTRMM']
                maxCEprecip = np.max(eachNode['cloudElementLatLonTRMM'][np.nonzero(eachNode['cloudElementLatLonTRMM'])])
                minCEprecip = np.min(eachNode['cloudElementLatLonTRMM'][np.nonzero(eachNode['cloudElementLatLonTRMM'])])
                mcsPrecip.append((eachNode['uniqueID'], minCEprecip, maxCEprecip))
        else:
            for eachMCC in finalMCCList:
                #get the info from the node
                for node in eachMCC:
                    eachNode = mccSearch.this_dict(node)
                    #find min and max precip
                    maxCEprecip = np.max(eachNode['cloudElementLatLonTRMM'][np.nonzero(eachNode['cloudElementLatLonTRMM'])])
                    minCEprecip = np.min(eachNode['cloudElementLatLonTRMM'][np.nonzero(eachNode['cloudElementLatLonTRMM'])])
                    mcsPrecip.append((eachNode['uniqueID'], minCEprecip, maxCEprecip))
                allmcsPrecip.append(mcsPrecip)
                mcsPrecip = []

    return mcsPrecip
#**********************************************************************************************************************
def precip_totals(finalMCCList):
    '''
    Purpose:: Precipitation totals associated with a cloud element

    Inputs:: finalMCCList: a list of dictionaries representing a list of nodes per a MCC

    Returns:: precipTotal: a floating-point number representing the total amount of precipitation associated
       with the feature

    '''
    precipTotal = 0.0
    cloudElementPrecip = 0.0
    mcsPrecip = []
    allmcsPrecip = []
    count = 0

    if finalMCCList:
        #print 'len finalMCCList is: ', len(finalMCCList)
        for eachMCC in finalMCCList:
            #get the info from the node
            for node in eachMCC:
                eachNode = mccSearch.this_dict(node)
                count += 1
                if count == 1:
                    prevHr = int(str(eachNode['cloudElementTime']).replace(' ', '')[-8:-6])

                currHr = int(str(eachNode['cloudElementTime']).replace(' ', '')[-8:-6])
                if prevHr == currHr:
                    cloudElementPrecip += eachNode['cloudElementPrecipTotal']
                else:
                    mcsPrecip.append((prevHr, cloudElementPrecip))
                    cloudElementPrecip = eachNode['cloudElementPrecipTotal']
                #last value in for loop
                if count == len(eachMCC):
                    mcsPrecip.append((currHr, cloudElementPrecip))

                precipTotal += eachNode['cloudElementPrecipTotal']
                prevHr = currHr

            mcsPrecip.append(('0', precipTotal))

            allmcsPrecip.append(mcsPrecip)
            precipTotal = 0.0
            cloudElementPrecip = 0.0
            mcsPrecip = []
            count = 0

        print 'allmcsPrecip ', allmcsPrecip

    return allmcsPrecip
#**********************************************************************************************************************
def shortest_duration(allMCCtimes):
    '''
    Purpose:: To determine the shortest MCC for the period

    Inputs:: allMCCtimes {MCCtimes, starttime, endtime, duration): a list of dictionaries of MCC temporal details for
        each MCC in the period considered

    Returns:: a floating-point number representing the minimum duration MCC detected

    Assumptions::

    '''

    return min([MCC['duration'] for MCC in allMCCtimes])
#**********************************************************************************************************************
def temporal_and_area_info_metric(finalMCCList):
    '''
    Purpose:: To provide information regarding the temporal properties of the MCCs found

    Inputs:: finalMCCList: a list of dictionaries representing a list of nodes representing a MCC

    Returns:: allMCCtimes: a list of dictionaries {mccTimes, starttime, endtime, duration, area} representing the MCC
       temporal details for each MCC in the period considered

    Assumptions::
        the final time hour --> the event lasted throughout that hr, therefore +1 to endtime
    '''
    #TODO: in real data edit this to use datetime
    mccTimes = []
    allMCCtimes = []
    mcsArea = []

    if finalMCCList:
        for eachMCC in finalMCCList:
            #get the info from the node
            for eachNode in eachMCC:
                mccTimes.append(mccSearch.this_dict(eachNode)['cloudElementTime'])
                mcsArea.append(mccSearch.this_dict(eachNode)['cloudElementArea'])

            #sort and remove duplicates
            mccTimes = list(set(mccTimes))
            mccTimes.sort()
            tdelta = mccTimes[1] - mccTimes[0]
            starttime = mccTimes[0]
            endtime = mccTimes[-1]
            duration = (endtime - starttime) + tdelta
            print 'starttime ', starttime, 'endtime ', endtime, 'tdelta ', tdelta, 'duration ', duration, 'MCSAreas ', mcsArea
            allMCCtimes.append({'MCCtimes':mccTimes, 'starttime':starttime, 'endtime':endtime, \
                'duration':duration, 'MCSArea': mcsArea})
            mccTimes = []
            mcsArea = []
    else:
        allMCCtimes = []
        tdelta = 0

    return allMCCtimes, tdelta
#**********************************************************************************************************************
