'''
# Wizard for running the mccSearch program
'''

import subprocess
import sys

import networkx as nx

# mccSearch modules
import mccSearch
import utils
import plotting
import metrics
import iomethods
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
    
    # For GrADs
    subprocess.call('export DISPLAY=:0.0', shell=True)

    # For first time working with the raw MERG zipped files
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

    # -------------------------------------------------------------------------------------------------
    # Getting started. Make it so number one!
    print ("-"*80)
    print "\t\t Starting the MCCSearch Analysis "
    print ("-"*80)
    print "\n -------------- Reading MERG Data ----------"
    mergImgs, timeList, LAT, LON, userVariables = iomethods.read_data('ch4', 'latitude', 'longitude', userVariables)
    print "\n -------------- findCloudElements ----------"
    CEGraph, _ = mccSearch.find_cloud_elements(mergImgs, timeList, LAT, LON, userVariables, graphVariables, userVariables.DIRS['TRMMdirName'])
    # theList = CEGraph.successors(node)
    # if the TRMMdirName wasnt entered for whatever reason, you can still get the TRMM data this way
    # CEGraph = mccSearch.findCloudElements(mergImgs,timeList)
    # allCETRMMList=mccSearch.findPrecipRate(DIRS['TRMMdirName'],timeList)
    # ----------------------------------------------------------------------------------------------
    print "\n -------------- findCloudClusters ----------"
    prunedGraph = mccSearch.find_cloud_clusters(CEGraph, userVariables, graphVariables)
    print "\n -------------- findMCCs ----------"
    MCCList, MCSList = mccSearch.find_MCC(prunedGraph, userVariables, graphVariables)
    # Now ready to perform various calculations/metrics
    print ("-"*80)
    print "\n -------------- METRICS ----------"
    print ("-"*80)
    # Some calculations/metrics that work that work
    print "creating the MCC userfile ", metrics.create_text_file(MCCList, 1, userVariables, graphVariables)
    print "creating the MCS userfile ", metrics.create_text_file(MCSList, 2, userVariables, graphVariables)
    plot_menu(MCCList, MCSList, userVariables.DIRS)

    #Let's get outta here! Engage!
    print ("-"*80)
# *********************************************************************************************************************


def plot_menu(MCCList, MCSList, DIRS):
    '''
    Purpose:: The flow of plots for the user to choose

    Input:: MCCList: a list of directories representing a list of nodes in the MCC
            MCSList: a list of directories representing a list of nodes in the MCS
            DIRS: a dictionary indicating the paths to data, both original and generated

    Output:: None
    '''
    option = display_plot_menu()
    while option != 0:
        #try:
        if option == 1:
            print "Generating Accumulated Rainfall from TRMM for the entire period ...\n"
            plotting.plot_accu_TRMM(MCSList, DIRS['mainDirStr'])
        if option == 2:
            startDateTime = raw_input("> Please enter the start date and time yyyy-mm-dd_hr:mm:ss format: \n")
            endDateTime = raw_input("> Please enter the end date and time yyyy-mm-dd_hr:mm:ss format: \n")
            print "Generating acccumulated rainfall between ", startDateTime, " and ", endDateTime, " ... \n"
            plotting.plot_accu_in_time_range(startDateTime, endDateTime, DIRS['mainDirStr'], 1.0)
        if option == 3:
            print "Generating area distribution plot ... \n"
            plotting.display_size(MCSList, DIRS['mainDirStr'])
        if option == 4:
            print "Generating precipitation and area distribution plot ... \n"
            plotting.display_precip(MCSList, DIRS['mainDirStr'])
        if option == 5:
            try:
                print "Generating histogram of precipitation for each time ... \n"
                plotting.plot_precip_histograms(MCSList, DIRS['mainDirStr'])
            except:
                pass
        # except:
        #     print "Invalid option. Please try again, enter 0 to exit \n"

        option = display_plot_menu()
    return
# *********************************************************************************************************************


def display_plot_menu():
    '''
    Purpose:: Display the plot Menu Options

    Input:: None

    Output:: option: an integer representing the choice of the user
    '''
    print "**************** PLOTS ************** \n"
    print "0. Exit \n"
    print "1. Accumulated TRMM precipitation \n"
    print "2. Accumulated TRMM precipitation between dates \n"
    print "3. Area distribution of the system over time \n"
    print "4. Precipitation and area distribution of the system \n"
    print "5. Histogram distribution of the rainfall in the area \n"
    option = int(raw_input("> Please enter your option for plots: \n"))
    return option
# *********************************************************************************************************************


def display_postprocessing_plot_menu():
    '''
    Purpose:: Display the plot Menu Options

    Input:: None

    Output:: option: an integer representing the choice of the user
    '''
    print "**************** POST PROCESSING PLOTS ************** \n"
    print "0. Exit \n"
    print "1. Map plots of the original MERG data \n"
    print "2. Map plots of the cloud elements using IR data \n"
    print "3. Map plots of the cloud elements rainfall accumulations using TRMM data \n"
    #print "4. Accumulated TRMM precipitation \n"
    #print "5. Accumulated TRMM precipitation between dates \n"

    option = int(raw_input("> Please enter your option for plots: \n"))
    return option
# *********************************************************************************************************************


def postprocessing_plot_menu(DIRS):
    '''
    Purpose:: The flow of plots for the user to choose

    Input:: DIRS a dictionary of directories
    #       DIRS={
    #          mainDirStr= "/directory/to/where/to/store/outputs"
    #          TRMMdirName = "/directory/to/the/TRMM/netCDF/files"
    #          CEoriDirName = "/directory/to/the/MERG/netCDF/files"
    #         }

    Output:: None
    '''

    TRMMCEdirName = DIRS['mainDirStr']+'/TRMMnetcdfCEs'
    CEdirName = DIRS['mainDirStr']+'/MERGnetcdfCEs'

    option = display_postprocessing_plot_menu()
    while option != 0:
        try:
            if option == 1:
                print "Generating images from the original MERG dataset ... \n"
                utils.post_processing_netcdf(3, DIRS['CEoriDirName'])
            if option == 2:
                print "Generating images from the cloud elements using MERG IR data ... \n"
                utils.post_processing_netcdf(1, CEdirName)
            if option == 3:
                print "Generating precipitation accumulation images from the cloud elements using TRMM data ... \n"
                utils.post_processing_netcdf(2, TRMMCEdirName)
            # if option == 4:
            #     print "Generating Accumulated TRMM rainfall from cloud elements for each MCS ... \n"
            #     featureType = int(raw_input("> Please enter type of MCS MCC-1 or MCS-2: \n"))
            #     if featureType == 1:
            #         filename = DIRS['mainDirStr']+'/textFiles/MCCPostProcessing.txt'
            #         try:
            #             if os.path.isfile(filename):
            #             #read each line as a list
            #         plotting.plotAccTRMM(DIRS['mainDirStr'])
            # if option == 5:
            #     mccSearch.plotAccuInTimeRange()
        except:
            print "Invalid option, please try again"
        option = display_postprocessing_plot_menu()
    return
# *********************************************************************************************************************

main()

