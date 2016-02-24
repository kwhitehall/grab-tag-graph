import json
import os
import subprocess

import networkx as nx

import iomethods
import utils


class UserVariables(object):
    # These will be assigned as the user determines which values they would like
    # then the global var's above will be assigned the according to these values
    def __init__(self, useJSON):  # useJSON=False): #self):
        data = None
        # Check if baselineDataDir.zip is unzipped, if not unzip it
        if not os.path.exists("../datadir"):
            subprocess.call('cd ../ ; unzip baselineDataDir.zip', shell=True)
        
        if useJSON is True:
            try:
                with open('../config.json') as f:
                    data = json.load(f)
                    f.close()

                self.LATMIN = float(data['LATMIN'])
                self.LATMAX = float(data['LATMAX'])
                self.LONMIN = float(data['LONMIN'])
                self.LONMAX = float(data['LONMAX'])
                self.startDateTime = data['startDateTime']
                self.endDateTime = data['endDateTime']
                self.XRES = float(data['XRES'])
                self.YRES = float(data['YRES'])
                self.TRES = float(data['TRES'])
                self.T_BB_MAX = float(data['T_BB_MAX'])
                self.T_BB_MIN = float(data['T_BB_MIN'])
                self.CONVECTIVE_FRACTION = float(data['CONVECTIVE_FRACTION'])
                self.MIN_MCS_DURATION = float(data['MIN_MCS_DURATION'])
                self.AREA_MIN = float(data['AREA_MIN'])
                self.MIN_OVERLAP = float(data['MIN_OVERLAP'])
                self.ECCENTRICITY_THRESHOLD_MAX = float(data['ECCENTRICITY_THRESHOLD_MAX'])
                self.ECCENTRICITY_THRESHOLD_MIN = float(data['ECCENTRICITY_THRESHOLD_MIN'])
                self.OUTER_CLOUD_SHIELD_AREA = float(data['OUTER_CLOUD_SHIELD_AREA'])
                self.INNER_CLOUD_SHIELD_AREA = float(data['INNER_CLOUD_SHIELD_AREA'])
                self.OUTER_CLOUD_SHIELD_TEMPERATURE = float(data['OUTER_CLOUD_SHIELD_TEMPERATURE'])
                self.INNER_CLOUD_SHIELD_TEMPERATURE = float(data['INNER_CLOUD_SHIELD_TEMPERATURE'])
                self.MINIMUM_DURATION = int(data['MINIMUM_DURATION'])
                self.MAXIMUM_DURATION = int(data['MAXIMUM_DURATION'])
                self.LAT_DISTANCE = float(data['LAT_DISTANCE'])    # The avg distance in km for 1deg lat for the region being considered
                self.LON_DISTANCE = float(data['LON_DISTANCE'])    # The avg distance in km for 1deg lon for the region being considered
                self.DIRS = data['DIRS']  

            except IOError, e:
                print "Config file not found! Using default variables..."

        else:
            self.LATMIN = '5.0'          # min latitude; -ve values in the SH e.g. 5S = -5
            self.LATMAX = '19.0'         # max latitude; -ve values in the SH e.g. 5S = -5 20.0
            self.LONMIN = '-5.0'         # min longitude; -ve values in the WH e.g. 59.8W = -59.8
            self.LONMAX = '9.0'          # min longitude; -ve values in the WH e.g. 59.8W = -59.8
            self.startDateTime = "200908310000"  # "yyyymmddhhss"
            self.endDateTime = "200908312100"
            self.XRES = 4.0              # x direction spatial resolution in km
            self.YRES = 4.0              # y direction spatial resolution in km
            self.TRES = 1                # temporal resolution in hrs
            self.LAT_DISTANCE = 111.0    # the avg distance in km for 1deg lat for the region being considered
            self.LON_DISTANCE = 111.0    # the avg distance in km for 1deg lon for the region being considered
            self.T_BB_MAX = 243          # warmest temp to allow (-30C to -55C according to Morel and Sensi 2002
            self.T_BB_MIN = 218          # cooler temp for the center of the system
            self.CONVECTIVE_FRACTION = 0.90  # the min temp/max temp that would be expected in a CE.. this is highly conservative (only a 10K difference)
            self.MIN_MCS_DURATION = 3    # minimum time for a MCS to exist
            self.AREA_MIN = 2400.0       # minimum area for CE criteria in km^2 according to Vila et al. (2008) is 2400
            self.MIN_OVERLAP = 10000.00  # km^2  see ref from Williams and Houze 1987 and Arnaud et al 1992
            self.ECCENTRICITY_THRESHOLD_MAX = 1.0       # tending to 1 is a circle e.g. hurricane,
            self.ECCENTRICITY_THRESHOLD_MIN = 0.70      # tending to 0 is a linear e.g. squall line
            self.OUTER_CLOUD_SHIELD_AREA = 80000.0      # km^2
            self.INNER_CLOUD_SHIELD_AREA = 30000.0      # km^2
            self.OUTER_CLOUD_SHIELD_TEMPERATURE = 233   # in K
            self.INNER_CLOUD_SHIELD_TEMPERATURE = 213   # in K
            self.MINIMUM_DURATION = 6    # min number of frames the MCC must exist for (assuming hrly frames, MCCs is 6hrs)
            self.MAXIMUM_DURATION = 24   # max number of frames the MCC can last for
            self.DIRS = {'mainDirStr': "../firstattempt",
                         'TRMMdirName': "../datadir/TRMM",
                         'CEoriDirName': "../datadir/MERG"}
            self.filelist = None

        self.STRUCTURING_ELEMENT = [[0, 1, 0],  # The matrix for determining the pattern for the contiguous boxes and must
                                    [1, 1, 1],  # have same rank of the matrix it is being compared against
                                    [0, 1, 0]   # criteria for determining cloud elements and edges
                                    ]

        self.check_lats(self.LATMIN, self.LATMAX)
        self.check_lons(self.LONMIN, self.LONMAX)
        self.check_times(self.startDateTime, self.endDateTime)
        self.check_dirs(self.DIRS['CEoriDirName'], self.DIRS['TRMMdirName'])
        self.ir_inputs(self.DIRS['CEoriDirName'], self.startDateTime, self.endDateTime)
        self.trmm_inputs(self.DIRS['TRMMdirName'], self.startDateTime, self.endDateTime)
        self.setup_all()
        
    def check_lats(self, LATMIN, LATMAX):
        '''
        Purpose:: check latitude limits
        '''

        if (float(self.LATMIN) < -90.):
            print "Bad latmin input! Check the config file. Now using default variables..."
            return False

        if (float(self.LATMAX) > 90.):
            print "Bad latmax input! Check the config file. Now using default variables..."
            return False   

        if (float(self.LATMIN) > float(self.LATMAX)):
            print 'LATMIN >= LATMAX; Reversing the latitude', self.LATMIN
            return False
        
        return True

    def check_lons(self, LONMIN, LONMAX):
        '''
        Purpose:: check longitude limits
        '''

        if (float(self.LONMIN) < -180.):
            print "Bad lonmin input! Check the config file. Now using default variables..."
            return False

        if (float(self.LONMAX) > 180.):
            print "Bad lonmax input! Check the config file. Now using default variables..."
            return False

        if (float(self.LONMIN) > float(self.LONMAX)):
            print 'LONMIN >= LONMAX; Reversing the longitude'
            return False

        return True

    def check_dirs(self, CEoriDirName, TRMMdirName):
        '''
        Purpose:: check the dirs
        '''

        if not os.path.exists(CEoriDirName):
            print "Error! MERG invalid path!"
            self.CEoriDirName = raw_input("> Please enter the directory to the MERG netCDF files: \n")
            return False
        if not TRMMdirName == "None":
            if not os.path.exists(TRMMdirName):
                print "Error: TRMM invalid path!"
                self.TRMMdirName = raw_input("> Please enter the location to the raw TRMM netCDF files: \n")
                return False
        return True

    def check_times(self, startDateTime, endDateTime):
        '''
        Purpose:: check the times
        '''

        # Check validity of time
        while utils.valid_date(startDateTime) != True:
            print "Invalid time entered for startDateTime!"
        while utils.valid_date(endDateTime) != True:
            print "Invalid time entered for endDateTime!"

        if (float(startDateTime) > float(endDateTime)):
            print "Bad start and end times input! Check the config file. Now using default variables..."
            return False
        return True

    def trmm_inputs(self, TRMMdirName, startDateTime, endDateTime):
        # Checks that inputs are ok
        if not TRMMdirName == "None":
            # Check if all the files exists in the TRMM directory entered
            test, _ = iomethods.check_for_files(TRMMdirName, startDateTime, endDateTime, 3, 'hour')
            if test is False:
                print "Error with files in the TRMM directory entered. Please check your files before restarting. "
                return False
    
        return True
    
    def ir_inputs(self, CEoriDirName, startDateTime, endDateTime):
        try:
            test, self.filelist = iomethods.check_for_files(CEoriDirName, startDateTime, endDateTime, 1, 'hour')
            if test is False:
                print "Error with files in the MERG directory entered. Please check your files before restarting. "
                return False
        except:
            print "..." 
        return True

    def setup_all(self):
        '''
        Purpose:: to configure the UserVariables

        Inputs::

        Outputs::

        Assumptions::
        '''

        if (self.XRES <= 0):
            print "Bad XRES input! Check the config file. Now using default variables..."
            return False
        
        if (self.YRES <= 0):
            print "Bad YRES input! Check the config file. Now using default variables..."
            return False
        
        if (self.TRES <= 0):
            print "Bad TRES input! Check the config file. Now using default variables..."
            return False

        if (self.LAT_DISTANCE > 111. or self.LAT_DISTANCE <= 0.):
            print "Bad LAT_DISTANCE input! Check the config file. Now using default variables..."
            return False
        
        if (self.LON_DISTANCE > 111. or self.LON_DISTANCE <= 0.):
            print "Bad LON_DISTANCE input! Check the config file. Now using default variables..."
            return False          

        if (self.T_BB_MIN > self.T_BB_MAX):
            print "Bad TBBMIN and MAX input! Check the config file. Switching max and min..."
            temp = self.T_BB_MIN
            self.T_BB_MIN = self.T_BB_MAX
            self.T_BB_MAX = temp

        if (self.CONVECTIVE_FRACTION < 0 or self.CONVECTIVE_FRACTION > 1):
            print "Bad convective fraction input. Now using default variables..."
            return False

        if (self.MIN_MCS_DURATION < 3):
            print "MIN_MCS_Duration too small. Minimum is 3. Now using default variables..."
            return False

        if (self.AREA_MIN <= 0):
            print "Bad AREA_MIN input. Now using default variables..."
            return False

        if (self.ECCENTRICITY_THRESHOLD_MAX < 0 or self.ECCENTRICITY_THRESHOLD_MAX > 1.):
            print "Bad ECCENTRICITY THRESHOLD MAX input. Check config file. Now using default variables..."
            return False
        
        if (self.ECCENTRICITY_THRESHOLD_MIN < 0 or self.ECCENTRICITY_THRESHOLD_MIN > 1.):
            print "Bad ECCENTRICITY THRESHOLD MIN input. Check config file. Now using default variables..."
            print False

        if (self.ECCENTRICITY_THRESHOLD_MIN > self.ECCENTRICITY_THRESHOLD_MAX):
            print "Bad ECCENTRICITY THRESHOLD MIN and MAX input! Check the config file. Now using default variables..."
            return False

        if (self.INNER_CLOUD_SHIELD_AREA > self.OUTER_CLOUD_SHIELD_AREA):
            print "Bad cloud shield areas input. Check config file. Now using default variables..."
            return False

        if (self.INNER_CLOUD_SHIELD_TEMPERATURE > self.OUTER_CLOUD_SHIELD_TEMPERATURE):
            print "Bad cloud shield temperatures. Check the config file. Now using default variables..."
            return False

        if (self.MINIMUM_DURATION > self.MAXIMUM_DURATION):
            print "Bad MIN and MAX DURATION inputs! Check the config file. Now using default variables..."
            return False

        return True


class GraphVariables(object):
    def __init__(self):
        self.edgeWeight = [1, 2, 3]                   # weights for the graph edges
        self.CLOUD_ELEMENT_GRAPH = nx.DiGraph()     # graph obj for the CEs meeting criteria
        self.PRUNED_GRAPH = nx.DiGraph()            # graph meeting the CC criteria
        

def define_graph_variables():
    '''
    Purpose:: to assign the data for the GraphVariables

    Inputs:: None

    Outputs:: a Graphvariables instance 

    Assumptions::
    '''

    graphVars = GraphVariables()
    return graphVars


