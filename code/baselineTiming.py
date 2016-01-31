'''
# running the program
'''
import unittest
import time
import sys
import subprocess

import networkx as nx

import iomethods
import mccSearch
import utils
import variables

def main():
    sys.setrecursionlimit(5000)

    CEGraph = nx.DiGraph()
    prunedGraph = nx.DiGraph()
    MCCList = []
    MCSList = []
    MCSMCCNodesList = []
    allMCSsList = []
    allCETRMMList = []
    DIRS = {}

    #for GrADs
    subprocess.call('export DISPLAY=:0.0', shell=True)

    #for first time working with the raw MERG zipped files
    # rawMERG = "/directory/to/the/raw/MERGfiles"
    # utils.preprocessing_merg(rawMERG)
    # ---------------------------------------------------------------------------------
    # ---------------------------------- user inputs --------------------------------------
    userVariables = variables.UserVariables(useJSON=False)
    graphVariables = variables.define_graph_variables()
    # ---------------------------------- end user inputs --------------------------------------
    
    # create main directory and file structure for storing intel
    userVariables.DIRS['mainDirStr'] = iomethods.create_main_directory(userVariables.DIRS['mainDirStr'])
    TRMMCEdirName = userVariables.DIRS['mainDirStr']+'/TRMMnetcdfCEs'
    CEdirName = userVariables.DIRS['mainDirStr']+'/MERGnetcdfCEs'

    unittestFile = open(userVariables.DIRS['mainDirStr']+'/textFiles/unittestResults.txt','wb')
    unittestFile.write("\n Timing results for "+userVariables.startDateTime+" to "+userVariables.endDateTime)

    #let's go!
    # time how long it takes to complete reading in the data
    print "\n Start the timer "
    starttime = time.time()
    print "\n -------------- Read MERG Data ----------"
    print "\n Start the timer for the data ingest process"
    readMergStart = time.time()
    mergImgs, timeList, LAT, LON, userVariables = iomethods.read_data('ch4','latitude','longitude', userVariables)
    readMergEnd = time.time()
    print "\n End the timer for the data ingest process"
    print "\n Total time to complete data ingest is %g seconds"%(readMergEnd - readMergStart)
    unittestFile.write("\n 1. Total time to complete data ingest is %g seconds"%(readMergEnd - readMergStart))
    print ("-"*80)

    print "\n -------------- TESTING findCloudElements ----------"
    print "\n Start the timer for findCloudElements process"
    findCEsStart = time.time()
    print "\n Using both MERG and TRMM simultaneously "
    #CEGraph = mccSearch.find_cloud_elements(mergImgs,timeList,userVariables.DIRS['mainDirStr'], LAT,LON,userVariables, graphVariables, userVariables.DIRS['TRMMdirName'])
    CEGraph,_ = mccSearch.find_cloud_elements(mergImgs,timeList, LAT, LON, userVariables, graphVariables, userVariables.DIRS['TRMMdirName'])
    findCEsEnd = time.time()
    print "\n Number of cloud elements found is: ", CEGraph.number_of_nodes()
    print "\n End the timer for findCloudElements process"
    print "\n Total time to complete finding cloud elements is %g seconds"%(findCEsEnd - findCEsStart)
    unittestFile.write("\n 2. Total time to complete finding cloud elements is %g seconds"%(findCEsEnd - findCEsStart))
    
    # #********* OR *******
    # # #timing each separately
    # CEGraph = mccSearch.find_cloud_elements(mergImgs,timeList,DIRS['mainDirStr'], LAT,LON)
    # findCEsEnd = time.time()
    # print "\n End the timer for findCloudElements process using MERG only"
    # print "\n Total time to complete finding cloud elements in MERG only is %g seconds"%(findCEsEnd - findCEsStart)
    # unittestFile.write("\n Total time to complete finding cloud elements in MERG only is %g seconds"%(findCEsEnd - findCEsStart))
    
    # # *** TRMM DATA SET TIMINGS *** 
    # print "\n Start the timer for findCloudElements process using TRMM only"
    # findCETRMMStart = time.time()
    # allCETRMMList = mccSearch.find_precip_rate(DIRS['TRMMdirName'],timeList)
    # findCETRMMEnd = time.time()
    # print "\n End the timer for findCloudElements process using TRMM only"
    # print "\n Total time to complete finding cloud elements in TRMM only is %g seconds"%(findCETRMMEnd - findCETRMMStart)
    # unittestFile.write("\n Total time to complete finding cloud elements in TRMM only is %g seconds"%(findCETRMMEnd - findCETRMMStart))
    # print "\n Number of cloud elements found is: ", CEGraph.number_of_nodes()
    # print "\n Total time to complete finding cloud elements is %g seconds"%(findCETRMMEnd - findCEsStart)
    # unittestFile.write("\n Total time to complete finding cloud elements is %g seconds"%(findCETRMMEnd - findCEsStart))
    # print ("-"*80)

    print "\n -------------- TESTING findCloudClusters ----------"
    print "\n Start the timer for findCloudClusters process"
    findCloudClustersStart = time.time()
    prunedGraph = mccSearch.find_cloud_clusters(CEGraph,userVariables,graphVariables)
    print "The number of nodes in the prunedGraph is: ", prunedGraph.number_of_nodes()
    findCloudClustersEnd = time.time()
    print "\n End the timer for the findCloudClusters process"
    print "\n Total time to complete finding the cloud clusters is %g seconds"%(findCloudClustersEnd - findCloudClustersStart)
    unittestFile.write("\n 3. Total time to complete finding the cloud clusters is %g seconds"%(findCloudClustersEnd - findCloudClustersStart))
    print "\n The CEGraph nodes are: %s "%CEGraph.nodes()
    unittestFile.write("\n The CEGraph nodes are: %s "%CEGraph.nodes())
    #GABE: Added for testing with checkMERGBaseline
    print "\n The pruned graph nodes are: %s "%prunedGraph.nodes()
    unittestFile.write("\n\n The pruned graph nodes are: %s "%prunedGraph.nodes())
     
    print ("-"*80)
    
    print "\n -------------- TESTING findMCCs ----------"
    print "\n Start the timer for the findMCCs process"
    findMCCStart = time.time()
    MCCList,MCSList = mccSearch.find_MCC(prunedGraph,userVariables,graphVariables)
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
    unittestFile.write("\n The CEGraph nodes are: %s" %CEGraph.nodes())
    unittestFile.write("\n The prunedGraph nodes are: %s" %prunedGraph.nodes())
    print ("-"*80)
    # TODO: report the domain

    unittestFile.close()

main()

