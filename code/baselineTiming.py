'''
# running the program
'''
import unittest
import time
import sys
import subprocess

import profile

import networkx as nx


import iomethods
import mccSearch
import utils


def main():
    sys.setrecursionlimit(5000)

    CEGraph = nx.DiGraph()
    prunedGraph = nx.DiGraph()
    MCCList =[]
    MCSList=[]
    MCSMCCNodesList =[]
    allMCSsList =[]
    allCETRMMList =[]
    DIRS={}

    #for GrADs
    subprocess.call('export DISPLAY=:0.0', shell=True)

    #for first time working with the raw MERG zipped files
    # rawMERG = "/directory/to/the/raw/MERGfiles"
    # utils.preprocessing_merg(rawMERG)
    # ---------------------------------------------------------------------------------
    # ---------------------------------- user inputs --------------------------------------
    DIRS['mainDirStr'] = "/home/student/grab-tag-graph/paperData/OUTPUTS"#"/directory/to/where/to/store/outputs"
    DIRS['TRMMdirName'] = "/home/student/grab-tag-graph/paperData/TRMM"#"/directory/to/the/TRMM/netCDF/files"
    DIRS['CEoriDirName'] = "/home/student/grab-tag-graph/paperData/MERG"#"/directory/to/the/MERG/netCDF/files"
    #get the dates for analysis
    startDateTime = "200609110000" #"yyyymmddhrmm"
    endDateTime = "200609121200"
    # ---------------------------------- end user inputs --------------------------------------
    # Checks that inputs are ok
    try:
        if not os.path.exists(DIRS['CEoriDirName']):
            print "Error! MERG invalid path!"
            DIRS['CEoriDirName'] = raw_input("> Please enter the directory to the MERG netCDF files: \n")
    except:
        print "..."

    try:
        if not os.path.exists(DIRS['TRMMdirName']):
            print "Error: TRMM invalid path!"
            DIRS['TRMMdirName'] = raw_input("> Please enter the location to the raw TRMM netCDF files: \n")
    except:
        pass

    try:
        if not os.path.exists(DIRS['CEoriDirName']):
            print "Error! MERG invalid path!"
            DIRS['CEoriDirName'] = raw_input("> Please enter the directory to the MERG netCDF files: \n")
    except:
        print "..."   

    #check validity of time
    while utils.valid_date(startDateTime) != True:
        print "Invalid time entered for startDateTime!"

    while utils.valid_date(endDateTime) != True:
        print "Invalid time entered for endDateTime!"
        
    #check if all the files exisits in the MERG and TRMM directories entered
    test,_ = iomethods.check_for_files(DIRS['TRMMdirName'], startDateTime, endDateTime, 3, 'hour')
    if test == False:
        print "Error with files in the TRMM directory entered. Please check your files before restarting. "
        return
    #test,filelist = iomethods.check_for_files(startDateTime, endDateTime, DIRS['CEoriDirName'],1)
    test,filelist = iomethods.check_for_files(DIRS['CEoriDirName'], startDateTime, endDateTime, 1, 'hour')

    if test == False:
        print "Error with files in the original MERG directory entered. Please check your files before restarting. "
        return
    # end checks 

    # create main directory and file structure for storing intel
    DIRS['mainDirStr'] = iomethods.create_main_directory(DIRS['mainDirStr'])
    TRMMCEdirName = DIRS['mainDirStr']+'/TRMMnetcdfCEs'
    CEdirName = DIRS['mainDirStr']+'/MERGnetcdfCEs'

    unittestFile = open(DIRS['mainDirStr']+'/unittestResults.txt','wb')
    unittestFile.write("\n Timing results for "+startDateTime+" to "+endDateTime)

    #let's go!
    # time how long it takes to complete reading in the data
    print "\n Start the timer "
    starttime = time.time()
    print "\n -------------- Read MERG Data ----------"
    print "\n Start the timer for the data ingest process"
    readMergStart = time.time()
    mergImgs, timeList, LAT, LON = iomethods.read_data(DIRS['CEoriDirName'],'ch4','latitude','longitude', filelist)
    readMergEnd = time.time()
    print "\n End the timer for the data ingest process"
    print "\n Total time to complete data ingest is %g seconds"%(readMergEnd - readMergStart)
    unittestFile.write("\n 1. Total time to complete data ingest is %g seconds"%(readMergEnd - readMergStart))
    print ("-"*80)

    print "\n -------------- TESTING findCloudElements ----------"
    print "\n Start the timer for findCloudElements process"
    findCEsStart = time.time()
    # ********* EITHER *********
    # print "\n Using both MERG and TRMM simultaneously "
    # CEGraph = mccSearch.find_cloud_elements(mergImgs,timeList,DIRS['mainDirStr'], LAT,LON,DIRS['TRMMdirName'])
    # findCEsEnd = time.time()
    # print "\n Number of cloud elements found is: ", CEGraph.number_of_nodes()
    # print "\n End the timer for findCloudElements process"
    # print "\n Total time to complete finding cloud elements is %g seconds"%(findCEsEnd - findCEsStart)
    # unittestFile.write("\n 2. Total time to complete finding cloud elements is %g seconds"%(findCEsEnd - findCEsStart))
    # print ("-"*80)
    # #********* OR *******
    # #timing each separately
    CEGraph = mccSearch.find_cloud_elements(mergImgs,timeList,DIRS['mainDirStr'], LAT,LON)
    findCEsEnd = time.time()
    print "\n End the timer for findCloudElements process using MERG only"
    print "\n Total time to complete finding cloud elements in MERG only is %g seconds"%(findCEsEnd - findCEsStart)
    unittestFile.write("\n Total time to complete finding cloud elements in MERG only is %g seconds"%(findCEsEnd - findCEsStart))
    print "\n Start the timer for findCloudElements process using TRMM only"
    findCETRMMStart = time.time()
    allCETRMMList = mccSearch.find_precip_rate(DIRS['TRMMdirName'],timeList)
    findCETRMMEnd = time.time()
    print "\n End the timer for findCloudElements process using TRMM only"
    print "\n Total time to complete finding cloud elements in TRMM only is %g seconds"%(findCEsEnd - findCEsStart)
    unittestFile.write("\n Total time to complete finding cloud elements in TRMM only is %g seconds"%(findCEsEnd - findCEsStart))
    print "\n Number of cloud elements found is: ", CEGraph.number_of_nodes()
    print "\n Total time to complete finding cloud elements is %g seconds"%(findCETRMMEnd - findCEsStart)
    unittestFile.write("\n Total time to complete finding cloud elements is %g seconds"%(findCETRMMEnd - findCEsStart))
    print ("-"*80)

    print "\n -------------- TESTING findCloudClusters ----------"
    print "\n Start the timer for findCloudClusters process"
    findCloudClustersStart = time.time()
    prunedGraph = mccSearch.find_cloud_clusters(CEGraph)
    print "The number of nodes in the prunedGraph is: ", prunedGraph.number_of_nodes()
    findCloudClustersEnd = time.time()
    print "\n End the timer for the findCloudClusters process"
    print "\n Total time to complete finding the cloud clusters is %g seconds"%(findCloudClustersEnd - findCloudClustersStart)
    unittestFile.write("\n 3. Total time to complete finding the cloud clusters is %g seconds"%(findCloudClustersEnd - findCloudClustersStart))
    print ("-"*80)
    
    print "\n -------------- TESTING findMCCs ----------"
    print "\n Start the timer for the findMCCs process"
    findMCCStart = time.time()
    MCCList,MCSList =  mccSearch.find_MCC(prunedGraph)
    #MCCList,MCSList = profile.runctx('mccSearch.find_MCC(prunedGraph)', globals(), {'mccSearch' : mccSearch, 'prunedGraph' : prunedGraph})
    print "\n MCC List has been acquired ", len(MCCList)
    print "\n MCS List has been acquired ", len(MCSList)
    findMCCEnd = time.time()
    print "\n End the timer for the findMCCs process"
    print "\n Total time to complete finding the MCCs is %g seconds"%(findMCCEnd - findMCCStart)
    unittestFile.write("\n 4. Total time to complete finding the MCCs is %g seconds"%(findMCCEnd - findMCCStart))
    print ("-"*80)
    #end the timer
    endtime = time.time()
    print ("*"*80)
    print "\n The entire evaluation took %g seconds to complete" %(endtime - starttime)
    unittestFile.write("\n The entire evaluation took %g seconds to complete" %(endtime - starttime))
    unittestFile.write("\n ----------------------------------------------------------------")
    unittestFile.write("\n Number of cloud elements found is: %d"%CEGraph.number_of_nodes())
    unittestFile.write("\n Number of edges (with the cloud elements) found is: %d"%CEGraph.number_of_edges())
    unittestFile.write("\n The number of nodes in the prunedGraph is: %d" %prunedGraph.number_of_nodes())
    unittestFile.write("\n The number of edges (with nodes) in the prunedGraph is: %d" %prunedGraph.number_of_edges())
    unittestFile.write("\n MCC List has been acquired %d" %len(MCCList))
    unittestFile.write("\n MCS List has been acquired %d" %len(MCSList))
    # TODO: report the domain

    unittestFile.close()

main()