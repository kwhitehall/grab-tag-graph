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

    #for GrADs
    subprocess.call('export DISPLAY=:0.0', shell=True)

    #for first time working with the raw MERG zipped files
    # rawMERG = "/directory/to/the/raw/MERGfiles"
    # utils.preprocessing_merg(rawMERG)
    # ---------------------------------------------------------------------------------
    # ---------------------------------- user inputs --------------------------------------
    DIRS['mainDirStr'] = "/Users/kwhitehall/Documents/capstones/usc2015/baselineTimings/test"#"/directory/to/where/to/store/outputs"
    DIRS['TRMMdirName'] = "/Users/kwhitehall/Documents/capstones/usc2015/baselineTimings/datadir/TRMM"#"/directory/to/the/TRMM/netCDF/files"
    DIRS['CEoriDirName'] = "/Users/kwhitehall/Documents/capstones/usc2015/baselineTimings/datadir/MERG"#"/directory/to/the/MERG/netCDF/files"
    #get the dates for analysis
    startDateTime = "200908310000" #"yyyymmddhrmm"
    endDateTime = "200908312100"
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
    
    # for doing some postprocessing with the clipped datasets instead of running the full program, e.g.
    # mccSearch.post_processing_netcdf(3,CEoriDirName)
    # mccSearch.post_processing_netcdf(2)
    # -------------------------------------------------------------------------------------------------

    #let's go!
    print "\n -------------- Read MERG Data ----------"
    mergImgs, timeList, LAT, LON = iomethods.read_data(DIRS['CEoriDirName'],'ch4','latitude','longitude', filelist)
    print ("-"*80)

    print "\n -------------- TESTING findCloudElements ----------"
    CEGraph = mccSearch.find_cloud_elements(mergImgs,timeList,DIRS['mainDirStr'], LAT,LON,DIRS['TRMMdirName'])
    #if the TRMMdirName wasnt entered for whatever reason, you can still get the TRMM data this way
    # CEGraph = mccSearch.findCloudElements(mergImgs,timeList,DIRS['mainDirStr'], LAT,LON)
    # allCETRMMList=mccSearch.findPrecipRate(DIRS['TRMMdirName'],timeList)
    # ----------------------------------------------------------------------------------------------
    print ("-"*80)
    print "number of nodes in CEGraph is: ", CEGraph.number_of_nodes()
    print ("-"*80)
    print "\n -------------- TESTING findCloudClusters ----------"
    prunedGraph = mccSearch.find_cloud_clusters(CEGraph)
    print ("-"*80)
    print "number of nodes in prunedGraph is: ", prunedGraph.number_of_nodes()
    print ("-"*80)
    #sys.exit()
    print "\n -------------- TESTING findMCCs ----------"
    MCCList,MCSList = mccSearch.find_MCC(prunedGraph)
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