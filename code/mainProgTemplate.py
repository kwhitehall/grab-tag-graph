'''
# running the program
'''

import networkx as nx
import subprocess

import iomethods
import mccSearch
import utils


def main():
    CEGraph = nx.DiGraph()
    prunedGraph = nx.DiGraph()
    MCCList =[]
    MCSList=[]
    MCSMCCNodesList =[]
    allMCSsList =[]
    allCETRMMList =[]
    DIRS={}

    # DIRS={
    #          mainDirStr= "/directory/to/where/to/store/outputs"
    #          TRMMdirName = "/directory/to/the/TRMM/netCDF/files"
    #          CEoriDirName = "/directory/to/the/MERG/netCDF/files"
    #         }

    #for GrADs
    subprocess.call('export DISPLAY=:0.0', shell=True)

    #for first time working with the raw MERG zipped files
    # rawMERG = "/directory/to/the/raw/MERGfiles"
    # utils.preprocessing_merg(rawMERG)
    # ---------------------------------------------------------------------------------


    #create main directory and file structure for storing intel
    DIRS['mainDirStr'] = utils.createMainDirectory(DIRS['mainDirStr'])
    TRMMCEdirName = DIRS['mainDirStr']+'/TRMMnetcdfCEs'
    CEdirName = DIRS['mainDirStr']+'/MERGnetcdfCEs'

    # for doing some postprocessing with the clipped datasets instead of running the full program, e.g.
    # mccSearch.post_processing_netcdf(3,CEoriDirName)
    # mccSearch.post_processing_netcdf(2)
    # -------------------------------------------------------------------------------------------------

    #let's go!
    print "\n -------------- Read MERG Data ----------"
    mergImgs, timeList, LAT, LON = iomethods.readMergData(DIRS['CEoriDirName'], filelist)
    print ("-"*80)

    print "\n -------------- TESTING findCloudElements ----------"
    CEGraph = mccSearch.findCloudElements(mergImgs,timeList,DIRS['mainDirStr'], LAT,LON,DIRS['TRMMdirName'])
    #if the TRMMdirName wasnt entered for whatever reason, you can still get the TRMM data this way
    # CEGraph = mccSearch.findCloudElements(mergImgs,timeList,DIRS['mainDirStr'], LAT,LON)
    # allCETRMMList=mccSearch.findPrecipRate(DIRS['TRMMdirName'],timeList)
    # ----------------------------------------------------------------------------------------------
    print ("-"*80)
    print "number of nodes in CEGraph is: ", CEGraph.number_of_nodes()
    print ("-"*80)
    print "\n -------------- TESTING findCloudClusters ----------"
    prunedGraph = mccSearch.findCloudClusters(CEGraph)
    print ("-"*80)
    print "number of nodes in prunedGraph is: ", prunedGraph.number_of_nodes()
    print ("-"*80)
    #sys.exit()
    print "\n -------------- TESTING findMCCs ----------"
    MCCList,MCSList = mccSearch.findMCC(prunedGraph)
    print ("-"*80)
    print "MCC List has been acquired ", len(MCCList)
    print "MCS List has been acquired ", len(MCSList)
    print ("-"*80)
    #now ready to perform various calculations/metrics
    print "\n -------------- TESTING METRICS ----------"

    #some calculations/metrics that work that work
    # print "creating the MCC userfile ", metrics.createTextFile(MCCList,1)
    # print "creating the MCS userfile ", metrics.createTextFile(MCSList,2)
    # MCCTimes, tdelta = metrics.temporalAndAreaInfoMetric(MCCList)
    # print "number of MCCs is: ", metrics.numberOfFeatures(MCCList)
    # print "longest duration is: ", metrics.longestDuration(MCCTimes), "hrs"
    # print "shortest duration is: ", metrics.shortestDuration(MCCTimes), "hrs"
    # print "Average duration is: ", metrics.averageDuration(MCCTimes), "hrs"
    # print "Average size is: ", mmetrics.averageFeatureSize(MCCList), "km^2"

    #some plots that work
    # plotting.plotAccTRMM(MCCList)
    # plotting.displayPrecip(MCCList)
    # plotting.plotAccuInTimeRange('yyyy-mm-dd_hh:mm:ss', 'yyyy-mm-dd_hh:mm:ss')
    # plotting.displaySize(MCCList)
    # plotting.displayPrecip(MCCList)
    # plotting.plotHistogram(MCCList)
    #
    print ("-"*80)

main()