import datetime
from datetime import timedelta, datetime

import numpy as np

import mccSearch

#****************************************************************** 
# 
#             METRICS FUNCTIONS FOR MERG.PY
# TODO: rewrite these metrics so that they use the data from the 
#   file instead of the graph as this reduce mem resources needed
#   
#
#******************************************************************
def averageDuration(allMCCtimes):
    '''
    Purpose:: To determine the average MCC length for the period

    Input:: list of dictionaries - allMCCtimes {MCCtimes, starttime, endtime, duration): a list of dictionaries
              of MCC temporal details for each MCC in the period considered

    Output::a floating-point representing the average duration of a MCC in the period
            
    Assumptions:: 

    '''

    return sum([MCC['duration'] for MCC in allMCCtimes], timedelta(seconds=0))/len(allMCCtimes)
#******************************************************************
def averageFeatureSize(finalMCCList): 
    '''
    Purpose:: To determine the average MCC size for the period

    Input:: a list of list of strings - finalMCCList: a list of list of nodes representing a MCC
    
    Output::a floating-point representing the average area of a MCC in the period
            
    Assumptions:: 

    '''
    thisMCC = 0.0
    thisMCCAvg = 0.0

    #for each node in the list, get the are information from the dictionary
    #in the graph and calculate the area
    for eachPath in finalMCCList:
        for eachNode in eachPath:
            thisMCC += mccSearch.thisDict(eachNode)['cloudElementArea']

        thisMCCAvg += (thisMCC/len(eachPath))
        thisMCC = 0.0

    #calcuate final average
    return thisMCCAvg/(len(finalMCCList))
#******************************************************************
def averageTime (allTimes):
    '''
    Purpose:: 
        To determine the average time in a list of datetimes 
        e.g. of use is finding avg starttime, 
    Input:: 
        allTimes: a list of datetimes representing all of a given event e.g. start time

    Output:: 
        a floating-point number representing the average of the times given

    '''
    avgTime = 0

    for aTime in allTimes:
        avgTime += aTime.second + 60*aTime.minute + 3600*aTime.hour

    if len(allTimes) > 1:
        avgTime /= len(allTimes)
    
    rez = str(avgTime/3600) + ' ' + str((avgTime%3600)/60) + ' ' + str(avgTime%60)
    return datetime.strptime(rez, "%H %M %S")
#******************************************************************
def commonFeatureSize(finalMCCList): 
    '''
    Purpose:: 
        To determine the common (mode) MCC size for the period

    Input:: 
        finalMCCList: a list of list of strings representing the list of nodes representing a MCC
    
    Output::
        a floating-point representing the average area of a MCC in the period
            
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
    return hist,bin_edges
#******************************************************************
def createTextFile(finalMCCList, identifier, MAINDIRECTORY, OUTER_CLOUD_SHIELD_AREA,TRES):
    '''
    Purpose:: 
        Create a text file with information about the MCS
        This function is expected to be especially of use regarding long term record checks

    Input:: 
        finalMCCList: a list of dictionaries representing a list of nodes representing a MCC
        identifier: an integer representing the type of list that has been entered...this is for creating file purposes
            1 - MCCList; 2- MCSList

    Output:: 
        a user readable text file with all information about each MCS
        a user readable text file with the summary of the MCS

    Assumptions:: 
    '''

    durations=0.0
    startTimes =[]
    endTimes=[]
    averagePropagationSpeed = 0.0
    speedCounter = 0
    maxArea =0.0
    amax = 0.0
    avgMaxArea =[]
    maxAreaCounter =0.0
    maxAreaTime=''
    eccentricity = 0.0
    firstTime = True
    matureFlag = True
    timeMCSMatures=''
    maxCEprecipRate = 0.0
    minCEprecipRate = 0.0
    averageArea = 0.0
    averageAreaCounter = 0
    durationOfMatureMCC = 0
    avgMaxPrecipRate = 0.0
    avgMaxPrecipRateCounter = 0
    avgMinPrecipRate = 0.0
    avgMinPrecipRateCounter = 0
    CEspeed = 0.0
    MCSspeed = 0.0
    MCSspeedCounter = 0
    MCSPrecipTotal = 0.0
    avgMCSPrecipTotalCounter = 0
    bigPtotal = 0.0
    bigPtotalCounter = 0
    allPropagationSpeeds =[]
    averageAreas =[]
    areaAvg = 0.0
    avgPrecipTotal = 0.0
    avgPrecipTotalCounter = 0
    avgMaxMCSPrecipRate = 0.0
    avgMaxMCSPrecipRateCounter = 0
    avgMinMCSPrecipRate = 0.0
    avgMinMCSPrecipRateCounter = 0
    minMax =[]
    avgPrecipArea = []
    location =[]
    avgPrecipAreaPercent = 0.0
    precipArea = 0.0
    precipAreaPercent = 0.0
    precipPercent =[]
    precipCounter = 0
    precipAreaAvg = 0.0
    minSpeed = 0.0
    maxSpeed =0.0

    if identifier == 1:
        MCSUserFile = open((MAINDIRECTORY+'/textFiles/MCCsUserFile.txt'),'wb')
        MCSSummaryFile = open((MAINDIRECTORY+'/textFiles/MCCSummary.txt'),'wb')
        MCSPostFile = open((MAINDIRECTORY+'/textFiles/MCCPostPrecessing.txt'),'wb')
    
    if identifier == 2:
        MCSUserFile = open((MAINDIRECTORY+'/textFiles/MCSsUserFile.txt'),'wb')
        MCSSummaryFile = open((MAINDIRECTORY+'/textFiles/MCSSummary.txt'),'wb')
        MCSPostFile = open((MAINDIRECTORY+'/textFiles/MCSPostPrecessing.txt'),'wb')

    for eachPath in finalMCCList:
        eachPath.sort(key=lambda nodeID:(len(nodeID.split('C')[0]), nodeID.split('C')[0], nodeID.split('CE')[1]))
        MCSPostFile.write("\n %s" %eachPath)

        startTime = mccSearch.thisDict(eachPath[0])['cloudElementTime']
        endTime = mccSearch.thisDict(eachPath[-1])['cloudElementTime']
        duration = (endTime - startTime) + timedelta(hours=TRES)
        
        # convert datatime duration to seconds and add to the total for the average duration of all MCS in finalMCCList
        durations += (duration.total_seconds()) 
        
        #durations += duration
        startTimes.append(startTime)
        endTimes.append(endTime)

        #get the precip info
        
        for eachNode in eachPath:

            thisNode = mccSearch.thisDict(eachNode)

            #set first time min "fake" values
            if firstTime == True:
                minCEprecipRate = thisNode['CETRMMmin']
                avgMinMCSPrecipRate += thisNode['CETRMMmin']
                firstTime = False

            # #calculate the speed
            # if thisNode['cloudElementArea'] >= OUTER_CLOUD_SHIELD_AREA:
            #   averagePropagationSpeed += findCESpeed(eachNode, eachPath)
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
                

            #calculations for only the mature stage
            if thisNode['nodeMCSIdentifier'] == 'M':
                #calculate average area of the maturity feature only 
                averageArea += thisNode['cloudElementArea']
                averageAreaCounter += 1
                durationOfMatureMCC +=1
                avgMaxPrecipRate += thisNode['CETRMMmax']
                avgMaxPrecipRateCounter += 1
                avgMinPrecipRate += thisNode['CETRMMmin']
                avgMinPrecipRateCounter += 1
                avgMaxMCSPrecipRate += thisNode['CETRMMmax']
                avgMaxMCSPrecipRateCounter += 1
                avgMinMCSPrecipRate += thisNode['CETRMMmin']
                avgMinMCSPrecipRateCounter += 1

                #the precip percentage (TRMM area/CE area)
                if thisNode['cloudElementArea'] >= 0.0 and thisNode['TRMMArea'] >= 0.0:
                    precipArea += thisNode['TRMMArea']
                    avgPrecipArea.append(thisNode['TRMMArea'])
                    avgPrecipAreaPercent += (thisNode['TRMMArea']/thisNode['cloudElementArea'])
                    precipPercent.append((thisNode['TRMMArea']/thisNode['cloudElementArea'])) 
                    precipCounter += 1

                # #system speed for only mature stage
                # CEspeed = findCESpeed(eachNode,eachPath)
                # if CEspeed > 0.0 :
                #   MCSspeed += CEspeed
                #   MCSspeedCounter += 1
                    
            #find accumulated precip
            if thisNode['cloudElementPrecipTotal'] > 0.0:
                MCSPrecipTotal += thisNode['cloudElementPrecipTotal']
                avgMCSPrecipTotalCounter +=1

        #A: calculate the average Area of the (mature) MCS
        if averageAreaCounter > 0: # and averageAreaCounter > 0:
            averageArea/= averageAreaCounter
            averageAreas.append(averageArea)

        #v: MCS speed 
        if MCSspeedCounter > 0: # and MCSspeed > 0.0:
            MCSspeed /= MCSspeedCounter
            
        #smallP_max: calculate the average max precip rate (mm/h)
        if avgMaxMCSPrecipRateCounter > 0 : #and avgMaxPrecipRate > 0.0:
            avgMaxMCSPrecipRate /= avgMaxMCSPrecipRateCounter
            
        #smallP_min: calculate the average min precip rate (mm/h)
        if avgMinMCSPrecipRateCounter > 0 : #and avgMinPrecipRate > 0.0:
            avgMinMCSPrecipRate /= avgMinMCSPrecipRateCounter
            
        #smallP_avg: calculate the average precipitation (mm hr-1)
        if MCSPrecipTotal > 0.0: # and avgMCSPrecipTotalCounter> 0:
            avgMCSPrecipTotal = MCSPrecipTotal/avgMCSPrecipTotalCounter
            avgPrecipTotal += avgMCSPrecipTotal
            avgPrecipTotalCounter += 1
            
        #smallP_total = MCSPrecipTotal
        #precip over the MCS lifetime prep for bigP_total
        if MCSPrecipTotal > 0.0: 
            bigPtotal += MCSPrecipTotal
            bigPtotalCounter += 1
            
        if maxArea > 0.0:
            avgMaxArea.append(maxArea)
            maxAreaCounter += 1

        #verage precipate area precentage (TRMM/CE area)
        if precipCounter > 0:
            avgPrecipAreaPercent /= precipCounter
            precipArea /= precipCounter


        #write stuff to file
        MCSUserFile.write("\n\n\nStarttime is: %s " %(str(startTime)))
        MCSUserFile.write("\nEndtime is: %s " %(str(endTime)))
        MCSUserFile.write("\nLife duration is %s hrs" %(str(duration)))
        MCSUserFile.write("\nTime of maturity is %s " %(timeMCSMatures))
        MCSUserFile.write("\nDuration mature stage is: %s " %durationOfMatureMCC*TRES)
        MCSUserFile.write("\nAverage area is: %.4f km^2 " %(averageArea))
        MCSUserFile.write("\nMax area is: %.4f km^2 " %(maxArea))
        MCSUserFile.write("\nMax area time is: %s " %(maxAreaTime))
        MCSUserFile.write("\nEccentricity at max area is: %.4f " %(eccentricity))
        MCSUserFile.write("\nCenter (lat,lon) at max area is: %.2f\t%.2f" %(location[0], location[1]))
        MCSUserFile.write("\nPropagation speed is %.4f " %(MCSspeed))
        MCSUserFile.write("\nMCS minimum preicip rate is %.4f mmh^-1" %(minCEprecipRate))
        MCSUserFile.write("\nMCS maximum preicip rate is %.4f mmh^-1" %(maxCEprecipRate))
        MCSUserFile.write("\nTotal precipitation during MCS is %.4f mm/lifetime" %(MCSPrecipTotal))
        MCSUserFile.write("\nAverage MCS precipitation is %.4f mm" %(avgMCSPrecipTotal))
        MCSUserFile.write("\nAverage MCS maximum precipitation is %.4f mmh^-1" %(avgMaxMCSPrecipRate))
        MCSUserFile.write("\nAverage MCS minimum precipitation is %.4f mmh^-1" %(avgMinMCSPrecipRate))
        MCSUserFile.write("\nAverage precipitation area is %.4f km^2 " %(precipArea))
        MCSUserFile.write("\nPrecipitation area percentage of mature system %.4f percent " %(avgPrecipAreaPercent*100))


        #append stuff to lists for the summary file
        if MCSspeed > 0.0:
            allPropagationSpeeds.append(MCSspeed)
            averagePropagationSpeed += MCSspeed
            speedCounter += 1

        #reset vars for next MCS in list
        aaveragePropagationSpeed = 0.0
        averageArea = 0.0
        averageAreaCounter = 0
        durationOfMatureMCC = 0
        MCSspeed = 0.0
        MCSspeedCounter = 0
        MCSPrecipTotal = 0.0
        avgMaxMCSPrecipRate =0.0
        avgMaxMCSPrecipRateCounter = 0
        avgMinMCSPrecipRate = 0.0
        avgMinMCSPrecipRateCounter = 0
        firstTime = True
        matureFlag = True
        avgMCSPrecipTotalCounter=0
        avgPrecipAreaPercent = 0.0
        precipArea = 0.0
        precipCounter = 0
        maxArea = 0.0
        maxAreaTime=''
        eccentricity = 0.0
        timeMCSMatures=''
        maxCEprecipRate = 0.0
        minCEprecipRate = 0.0
        location =[]

    #LD: average duration
    if len(finalMCCList) > 1:
        durations /= len(finalMCCList)
        durations /= 3600.0 #convert to hours
    
        #A: average area
        areaAvg = sum(averageAreas)/ len(finalMCCList)
    #create histogram plot here
    #if len(averageAreas) > 1:
    #   plotHistogram(averageAreas, "Average Area [km^2]", "Area [km^2]")

    #Amax: average maximum area
    if maxAreaCounter > 0.0: #and avgMaxArea > 0.0 : 
        amax = sum(avgMaxArea)/ maxAreaCounter
        #create histogram plot here
        #if len(avgMaxArea) > 1:
        #   plotHistogram(avgMaxArea, "Maximum Area [km^2]", "Area [km^2]")

    #v_avg: calculate the average propagation speed 
    if speedCounter > 0 :  # and averagePropagationSpeed > 0.0
        averagePropagationSpeed /= speedCounter
    
    #bigP_min: calculate the min rate in mature system
    if avgMinPrecipRate >  0.0: # and avgMinPrecipRateCounter > 0.0:
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
        #   plotHistogram(avgPrecipArea, "Average Rainfall Area [km^2]", "Area [km^2]")
        

    sTime = str(averageTime(startTimes))
    eTime = str(averageTime(endTimes))
    if len (allPropagationSpeeds) > 1:
        maxSpeed = max(allPropagationSpeeds)
        minSpeed = min(allPropagationSpeeds)
    
    #write stuff to the summary file
    MCSSummaryFile.write("\nNumber of features is %d " %(len(finalMCCList)))
    MCSSummaryFile.write("\nAverage duration is %.4f hrs " %(durations))
    MCSSummaryFile.write("\nAverage startTime is %s " %(sTime[-8:]))
    MCSSummaryFile.write("\nAverage endTime is %s " %(eTime[-8:]))
    MCSSummaryFile.write("\nAverage size is %.4f km^2 " %(areaAvg))
    MCSSummaryFile.write("\nAverage precipitation area is %.4f km^2 " %(precipAreaAvg))
    MCSSummaryFile.write("\nAverage maximum size is %.4f km^2 " %(amax))
    MCSSummaryFile.write("\nAverage propagation speed is %.4f ms^-1" %(averagePropagationSpeed))
    MCSSummaryFile.write("\nMaximum propagation speed is %.4f ms^-1 " %(maxSpeed))
    MCSSummaryFile.write("\nMinimum propagation speed is %.4f ms^-1 " %(minSpeed))
    MCSSummaryFile.write("\nAverage minimum precipitation rate is %.4f mmh^-1" %(avgMinPrecipRate))
    MCSSummaryFile.write("\nAverage maximum precipitation rate is %.4f mm h^-1" %(avgMaxPrecipRate))
    MCSSummaryFile.write("\nAverage precipitation is %.4f mm h^-1 " %(avgPrecipTotal))
    MCSSummaryFile.write("\nAverage total precipitation during MCSs is %.4f mm/LD " %(bigPtotal))
    MCSSummaryFile.write("\nAverage precipitation area percentage is %.4f percent " %(precipAreaPercent))


    MCSUserFile.close
    MCSSummaryFile.close
    MCSPostFile.close
    return
#******************************************************************
def findCESpeed(node, MCSList, theList):
    '''
    Purpose:: 
        To determine the speed of the CEs uses vector displacement delta_lat/delta_lon (y/x)

    Input:: 
        node: a string representing the CE
        MCSList: a list of strings representing the feature

    Output::
        CEspeed: a floating-point number representing the speed of the CE 

    '''

    delta_lon =0.0
    delta_lat =0.0
    CEspeed =[]
    theSpeed = 0.0
    
    nodeLatLon=mccSearch.thisDict(node)['cloudElementCenter']

    
    for aNode in theList:
        if aNode in MCSList:
            #if aNode is part of the MCSList then determine distance
            aNodeLatLon = mccSearch.thisDict(aNode)['cloudElementCenter']
            #calculate CE speed
            #checking the lats
            # nodeLatLon[0] += 90.0
            # aNodeLatLon[0] += 90.0
            # delta_lat = (nodeLatLon[0] - aNodeLatLon[0]) 
            delta_lat = ((mccSearch.thisDict(node)['cloudElementCenter'][0] +90.0) - (mccSearch.thisDict(aNode)['cloudElementCenter'][0]+90.0))
            # nodeLatLon[1] += 360.0
            # aNodeLatLon[1] += 360.0
            # delta_lon = (nodeLatLon[1] - aNodeLatLon[1]) 
            delta_lon = ((mccSearch.thisDict(node)['cloudElementCenter'][1]+360.0) - (mccSearch.thisDict(aNode)['cloudElementCenter'][1]+360.0))
            #failsafe for movement only in one dir
            if delta_lat == 0.0:
                delta_lat = 1.0

            if delta_lon == 0.0:
                delta_lon = 1.0
            
            try:
                theSpeed = abs((((delta_lat/delta_lon)*LAT_DISTANCE*1000)/(TRES*3600))) #convert to s --> m/s
            except:
                theSpeed = 0.0
            
            CEspeed.append(theSpeed)

            # print "~~~ ", thisDict(aNode)['uniqueID']
            # print "*** ", nodeLatLon, thisDict(node)['cloudElementCenter']
            # print "*** ", aNodeLatLon, thisDict(aNode)['cloudElementCenter']
            
    if not CEspeed:
        return 0.0
    else:
        return min(CEspeed) 
#******************************************************************
def longestDuration(allMCCtimes):
    '''
    Purpose:: 
        To determine the longest MCC for the period

    Input:: 
        allMCCtimes: a list of dictionaries {MCCtimes, starttime, endtime, duration, area} representing a list of dictionaries
            of MCC temporal details for each MCC in the period considered

    Output::
        an integer - lenMCC: representing the duration of the longest MCC found
           a list of strings - longestMCC: representing the nodes of longest MCC

    Assumptions:: 

    '''

    # MCCList = []
    # lenMCC = 0
    # longestMCC =[]

    # #remove duplicates
    # MCCList = list(set(finalMCCList))

    # longestMCC = max(MCCList, key = lambda tup:len(tup))
    # lenMCC = len(longestMCC)

    # return lenMCC, longestMCC

    return max([MCC['duration'] for MCC in allMCCtimes])
#******************************************************************
def numberOfFeatures(finalMCCList):
    '''
    Purpose:: 
        To count the number of MCCs found for the period

    Input:: 
        finalMCCList: a list of list of strings representing a list of list of nodes representing a MCC
    
    Output::
        an integer representing the number of MCCs found

    '''
    return len(finalMCCList)
#******************************************************************
def precipMaxMin(finalMCCList):
    '''
    TODO: this doesnt work the np.min/max function seems to be not working with the nonzero option..possibly a problem upstream with cloudElementLatLonTRMM
    Purpose:: 
        Precipitation maximum and min rates associated with each CE in MCS
    Input:: 
        finalMCCList: a list of dictionaries representing a list of nodes representing a MCC

    Output::
        MCSPrecip: a list indicating max and min rate for each CE identified

    '''
    maxCEprecip = 0.0
    minCEprecip =0.0
    MCSPrecip=[]
    allMCSPrecip =[]


    if finalMCCList:
        if type(finalMCCList[0]) is str: # len(finalMCCList) == 1:
            for node in finalMCCList:
                eachNode = mccSearch.thisDict(node)
                CETRMM = eachNode['cloudElementLatLonTRMM']

                # print "all ", np.min(CETRMM[np.nonzero(CETRMM)])
                # print "minCEprecip ", np.min(eachNode['cloudElementLatLonTRMM']) #[np.nonzero(eachNode['cloudElementLatLonTRMM'])])

                # print "maxCEprecip ", np.max(eachNode['cloudElementLatLonTRMM'][np.nonzero(eachNode['cloudElementLatLonTRMM'])])
                # sys.exit()
                maxCEprecip = np.max(eachNode['cloudElementLatLonTRMM'][np.nonzero(eachNode['cloudElementLatLonTRMM'])])
                minCEprecip = np.min(eachNode['cloudElementLatLonTRMM'][np.nonzero(eachNode['cloudElementLatLonTRMM'])])
                MCSPrecip.append((eachNode['uniqueID'],minCEprecip, maxCEprecip))
            
        else:
            for eachMCC in finalMCCList:
                #get the info from the node
                for node in eachMCC: 
                    eachNode=mccSearch.thisDict(node)
                    #find min and max precip
                    maxCEprecip =  np.max(eachNode['cloudElementLatLonTRMM'][np.nonzero(eachNode['cloudElementLatLonTRMM'])])
                    minCEprecip =  np.min(eachNode['cloudElementLatLonTRMM'][np.nonzero(eachNode['cloudElementLatLonTRMM'])])
                    MCSPrecip.append((eachNode['uniqueID'],minCEprecip, maxCEprecip))
                allMCSPrecip.append(MCSPrecip)
                MCSPrecip =[]
     
    return MCSPrecip
#******************************************************************
def precipTotals(finalMCCList):
    '''
    Purpose:: 
        Precipitation totals associated with a cloud element

    Input:: 
        finalMCCList: a list of dictionaries representing a list of nodes representing a MCC

    Output:: 
        precipTotal: a floating-point number representing the total amount of precipitation associated 
            with the feature
    '''
    precipTotal = 0.0
    CEprecip =0.0
    MCSPrecip=[]
    allMCSPrecip =[]
    count = 0

    if finalMCCList:
        #print "len finalMCCList is: ", len(finalMCCList)
        for eachMCC in finalMCCList:
            #get the info from the node
            for node in eachMCC:
                eachNode=mccSearch.thisDict(node)
                count += 1
                if count == 1:
                    prevHr = int(str(eachNode['cloudElementTime']).replace(" ", "")[-8:-6])
                
                currHr =int(str(eachNode['cloudElementTime']).replace(" ", "")[-8:-6])
                if prevHr == currHr:
                    CEprecip += eachNode['cloudElementPrecipTotal'] 
                else:
                    MCSPrecip.append((prevHr,CEprecip))
                    CEprecip = eachNode['cloudElementPrecipTotal'] 
                #last value in for loop
                if count == len(eachMCC):
                    MCSPrecip.append((currHr, CEprecip))

                precipTotal += eachNode['cloudElementPrecipTotal'] 
                prevHr = currHr

            MCSPrecip.append(('0',precipTotal))
            
            allMCSPrecip.append(MCSPrecip)
            precipTotal =0.0
            CEprecip = 0.0
            MCSPrecip = []
            count = 0

        print "allMCSPrecip ", allMCSPrecip

    return allMCSPrecip
#******************************************************************
def shortestDuration(allMCCtimes):
    '''
    Purpose:: To determine the shortest MCC for the period

    Input:: list of dictionaries - allMCCtimes {MCCtimes, starttime, endtime, duration): a list of dictionaries
              of MCC temporal details for each MCC in the period considered

    Output::an integer - lenMCC: representing the duration of the shortest MCC found
            a list of strings - longestMCC: representing the nodes of shortest MCC

    Assumptions:: 

    '''
    # lenMCC = 0
    # shortestMCC =[]
    # MCCList =[]
    
    # #remove duplicates
    # MCCList = list(set(finalMCCList))

    # shortestMCC = min(MCCList, key = lambda tup:len(tup))
    # lenMCC = len(shortestMCC)

    # return lenMCC, shortestMCC
    return min([MCC['duration'] for MCC in allMCCtimes])
#******************************************************************
def temporalAndAreaInfoMetric(finalMCCList):
    '''
    Purpose:: 
        To provide information regarding the temporal properties of the MCCs found

    Input:: 
        finalMCCList: a list of dictionaries representing a list of nodes representing a MCC
    
    Output:: 
        allMCCtimes: a list of dictionaries {MCCtimes, starttime, endtime, duration, area} representing a list of dictionaries
            of MCC temporal details for each MCC in the period considered

    Assumptions:: 
        the final time hour --> the event lasted throughout that hr, therefore +1 to endtime
    '''
    #TODO: in real data edit this to use datetime
    #starttime =0
    #endtime =0
    #duration = 0
    MCCtimes =[]
    allMCCtimes=[]
    MCSArea =[]
    
    if finalMCCList:
        for eachMCC in finalMCCList:
            #get the info from the node
            for eachNode in eachMCC:
                MCCtimes.append(mccSearch.thisDict(eachNode)['cloudElementTime'])
                MCSArea.append(mccSearch.thisDict(eachNode)['cloudElementArea'])
            
            #sort and remove duplicates 
            MCCtimes=list(set(MCCtimes))
            MCCtimes.sort()
            tdelta = MCCtimes[1] - MCCtimes[0]
            starttime = MCCtimes[0]
            endtime = MCCtimes[-1]
            duration = (endtime - starttime) + tdelta
            print "starttime ", starttime, "endtime ", endtime, "tdelta ", tdelta, "duration ", duration, "MCSAreas ", MCSArea
            allMCCtimes.append({'MCCtimes':MCCtimes, 'starttime':starttime, 'endtime':endtime, 'duration':duration, 'MCSArea': MCSArea})
            MCCtimes=[]
            MCSArea=[]
    else:
        allMCCtimes =[]
        tdelta = 0 

    return allMCCtimes, tdelta
#******************************************************************
