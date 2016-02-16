'''
# running the program
'''
import subprocess
import sys

import networkx as nx

import iomethods
import mccSearch
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
    DIRS={}

    # For GrADs
    subprocess.call('export DISPLAY=:0.0', shell=True)

    userVariables = variables.define_user_variables(useJSON=True)
    graphVariables = variables.define_graph_variables()

    # Create main directory and file structure for storing intel
    userVariables.DIRS['mainDirStr'] = iomethods.create_main_directory(userVariables.DIRS['mainDirStr'])
    TRMMCEdirName = userVariables.DIRS['mainDirStr']+'/TRMMnetcdfCEs'
    CEdirName = userVariables.DIRS['mainDirStr']+'/MERGnetcdfCEs'
    
    # Let's go!
    # ----------------------------------CORE GTG STEPS------------------------------------------------
    print "\n -------------- Read MERG Data ----------"
    mergImgs, timeList, LAT, LON, userVariables = iomethods.read_data('ch4', 'latitude', 'longitude', userVariables)
    print ("-"*80)
    print "\n -------------- TESTING findCloudElements ----------"
    CEGraph, _ = mccSearch.find_cloud_elements(mergImgs, timeList, LAT, LON, userVariables, graphVariables, userVariables.DIRS['TRMMdirName'])
    # #********* OR *******
    # CEGraph = mccSearch.find_cloud_elements(mergImgs,timeList,LAT, LON, userVariables, graphVariables)
    # allCETRMMList = mccSearch.find_precip_rate(userVariables.DIRS['TRMMdirName'],timeList)
    # #********************
    print ("-"*80)
    print "number of nodes in CEGraph is: ", CEGraph.number_of_nodes()
    print ("-"*80)
    print "\n -------------- TESTING findCloudClusters ----------"
    prunedGraph = mccSearch.find_cloud_clusters(CEGraph, userVariables, graphVariables)
    print ("-"*80)
    print "number of nodes in prunedGraph is: ", prunedGraph.number_of_nodes()
    print ("-"*80)
    print "\n -------------- TESTING findMCCs ----------"
    MCCList, MCSList = mccSearch.find_MCC(prunedGraph, userVariables, graphVariables)
    print ("-"*80)
    print "MCC List has been acquired ", len(MCCList)
    print "MCS List has been acquired ", len(MCSList)
    print ("-"*80)
    # ---------------------------------END CORE GTG STEPS----------------------------------------------
    # Now ready to perform various calculations/metrics
    print "\n -------------- TESTING METRICS ----------"
    # Some calculations/metrics that work that work
    # print "creating the MCC userfile ", metrics.createTextFile(MCCList,1)
    # print "creating the MCS userfile ", metrics.createTextFile(MCSList,2)
    # MCCTimes, tdelta = metrics.temporalAndAreaInfoMetric(MCCList)
    # print "number of MCCs is: ", metrics.numberOfFeatures(MCCList)
    # print "longest duration is: ", metrics.longestDuration(MCCTimes), "hrs"
    # print "shortest duration is: ", metrics.shortestDuration(MCCTimes), "hrs"
    # print "Average duration is: ", metrics.averageDuration(MCCTimes), "hrs"
    # print "Average size is: ", metrics.averageFeatureSize(MCCList), "km^2"

    # Some plots that work
    # plotting.plotAccTRMM(MCCList)
    # plotting.displayPrecip(MCCList)
    # plotting.plotAccuInTimeRange('yyyy-mm-dd_hh:mm:ss', 'yyyy-mm-dd_hh:mm:ss')
    # plotting.displaySize(MCCList)
    # plotting.displayPrecip(MCCList)
    # plotting.plotHistogram(MCCList)
    #
    print ("-"*80)

main()
