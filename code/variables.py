import os
import networkx as nx
import json
import utils
import iomethods

class UserVariables(object):
    # these will be assigned as the user determines which values they would like 
    # then the global var's above will be assigned the according to these values
    def __init__(self):
        self.LATMIN = '5.0'          #min latitude; -ve values in the SH e.g. 5S = -5
        self.LATMAX = '19.0'         #max latitude; -ve values in the SH e.g. 5S = -5 20.0
        self.LONMIN = '-5.0'         #min longitude; -ve values in the WH e.g. 59.8W = -59.8 
        self.LONMAX = '9.0'          #min longitude; -ve values in the WH e.g. 59.8W = -59.8  
        self.XRES = 4.0              #x direction spatial resolution in km
        self.YRES = 4.0              #y direction spatial resolution in km
        self.TRES = 1                #temporal resolution in hrs
        self.LAT_DISTANCE = 111.0    #the avg distance in km for 1deg lat for the region being considered                     
        self.LON_DISTANCE = 111.0    #the avg distance in km for 1deg lon for the region being considered                     
        self.STRUCTURING_ELEMENT = [[0, 1, 0],
                                    [1, 1, 1],
                                    [0, 1, 0]
                                    ] #the matrix for determining the pattern for the contiguous boxes and must                        
#have same rank of the matrix it is being compared against
#criteria for determining cloud elements and edges                                  
        self.T_BB_MAX = 243          #warmest temp to allow (-30C to -55C according to Morel and Sensi 2002
        self.T_BB_MIN = 218          #cooler temp for the center of the system
        self.CONVECTIVE_FRACTION = 0.90 #the min temp/max temp that would be expected in a CE.. this is highly conservative (only a 10K difference)
        self.MIN_MCS_DURATION = 3    #minimum time for a MCS to exist
        self.AREA_MIN = 2400.0       #minimum area for CE criteria in km^2 according to Vila et al. (2008) is 2400
        self.MIN_OVERLAP = 10000.00  #km^2  see ref from Williams and Houze 1987 and Arnaud et al 1992                            
        self.ECCENTRICITY_THRESHOLD_MAX = 1.0       #tending to 1 is a circle e.g. hurricane,
        self.ECCENTRICITY_THRESHOLD_MIN = 0.70      #tending to 0 is a linear e.g. squall line
        self.OUTER_CLOUD_SHIELD_AREA = 80000.0      #km^2
        self.INNER_CLOUD_SHIELD_AREA = 30000.0      #km^2
        self.OUTER_CLOUD_SHIELD_TEMPERATURE = 233   #in K
        self.INNER_CLOUD_SHIELD_TEMPERATURE = 213   #in K                                                                 
        self.MINIMUM_DURATION = 6    #min number of frames the MCC must exist for (assuming hrly frames, MCCs is 6hrs) 
        self.MAXIMUM_DURATION = 24   #max number of framce the MCC can last for   
        self.DIRS = {'mainDirStr': "/Users/kwhitehall/Documents/capstones/usc2015/haig", \
                     'TRMMdirName':"/Users/kwhitehall/Documents/capstones/usc2015/paperData/TRMM", \
                     'CEoriDirName': "/Users/kwhitehall/Documents/capstones/usc2015/paperData/MERG"}
        self.startDateTime = "200609110000"
        self.endDateTime = "200609121200"
        self.filelist = None


    def setupJSON(self):
        '''
        Purpose:: to configure the UserVariables

        Inputs::

        Outputs::

        Assumptions::
        '''
        data = None
        try:
            with open('../config.json') as f:
                data = json.load(f)
                f.close()
        except IOError, e:
            print "Config file not found! Using default variables..."
            return False
        try:
            self.LATMIN = float(data['LATMIN'])
            if (self.LATMIN < -90. or self.LATMIN > 90.):
                print "Bad latmin input! Check the config file. Now using default variables..."
                return False

            self.LATMAX = float(data['LATMAX'])
            if (self.LATMAX < -90. or self.LATMAX > 90. or self.LATMAX < self.LATMIN):
                print "Bad latmax input! Check the config file. Now using default variables..."
                return False

            self.LONMIN = float(data['LONMIN'])
            if (self.LONMIN < -180. or self.LONMIN > 180.):
                print "Bad lonmin input! Check the config file. Now using default variables..."
                return False
            
            self.LONMAX = float(data['LONMAX'])
            if (self.LONMAX < -180. or self.LONMAX > 180. or self.LONMAX < self.LONMIN):
                print "Bad lonmax input! Check the config file. Now using default variables..."
                return False
            
            self.XRES = float(data['XRES'])
            if (self.XRES <= 0):
                print "Bad XRES input! Check the config file. Now using default variables..."
                return False
            
            self.YRES = float(data['YRES'])
            if (self.YRES <= 0):
                print "Bad YRES input! Check the config file. Now using default variables..."
                return False
            
            self.TRES = float(data['TRES'])
            if (self.TRES <= 0):
                print "Bad TRES input! Check the config file. Now using default variables..."
                return False
            
            self.LAT_DISTANCE = float(data['LAT_DISTANCE'])
            # if (self.LAT_DISTANCE > 111.):
            #     print "Bad LAT_DISTANCE input! Check the config file. Now using default variables..."
            #     return False
            
            self.LON_DISTANCE = float(data['LON_DISTANCE'])
            # if (self.LON_DISTANCE > 111.):
            #     print "Bad LON_DISTANCE input! Check the config file. Now using default variables..."
            #     return False          

            self.STRUCTURING_ELEMENT = data['STRUCTURING_ELEMENT'] 
            self.T_BB_MAX = float(data['T_BB_MAX'])
            self.T_BB_MIN = float(data['T_BB_MIN'])
            if (self.T_BB_MIN > self.T_BB_MAX):
                print "Bad TBBMIN and MAX input! Check the config file. Switching max and min..."
                temp = self.T_BB_MIN
                self.T_BB_MIN = self.T_BB_MAX
                self.T_BB_MAX = temp

            self.CONVECTIVE_FRACTION = float(data['CONVECTIVE_FRACTION'])
            if (self.CONVECTIVE_FRACTION < 0 or self.CONVECTIVE_FRACTION > 1):
                print "Bad convective fraction input. Now using default variables..."
                return False

            self.MIN_MCS_DURATION = float(data['MIN_MCS_DURATION'])
            if (self.MIN_MCS_DURATION < 3):
                print "MIN_MCS_Duration too small. Minimum is 3. Now using default variables..."
                return False

            self.AREA_MIN = float(data['AREA_MIN'])
            if (self.AREA_MIN <= 0):
        		print "Bad AREA_MIN input. Now using default variables..."
        		return False

            self.MIN_OVERLAP = float(data['MIN_OVERLAP'])
            self.ECCENTRICITY_THRESHOLD_MAX = float(data['ECCENTRICITY_THRESHOLD_MAX'])
            self.ECCENTRICITY_THRESHOLD_MIN = float(data['ECCENTRICITY_THRESHOLD_MIN'])
            if (self.ECCENTRICITY_THRESHOLD_MAX < 0 or self.ECCENTRICITY_THRESHOLD_MAX > 1):
                print "Bad ECCENTRICITY THRESHOLD MAX input. Check config file. Now using default variables..."
                return False
            if (self.ECCENTRICITY_THRESHOLD_MIN < 0 or self.ECCENTRICITY_THRESHOLD_MIN > 1):
        		print "Bad ECCENTRICITY THRESHOLD MIN input. Check config file. Now using default variables..."
        		print False

            if (self.ECCENTRICITY_THRESHOLD_MIN > self.ECCENTRICITY_THRESHOLD_MAX):
                print "Bad ECCENTRICITY THRESHOLD MIN and MAX input! Check the config file. Now using default variables..."
                return False;        

            self.OUTER_CLOUD_SHIELD_AREA = float(data['OUTER_CLOUD_SHIELD_AREA'])
            self.INNER_CLOUD_SHIELD_AREA = float(data['INNER_CLOUD_SHIELD_AREA'])
            if (self.INNER_CLOUD_SHIELD_AREA > self.OUTER_CLOUD_SHIELD_AREA):
        		print "Bad cloud shield areas input. Check config file. Now using default variables..."
        		return False

            self.OUTER_CLOUD_SHIELD_TEMPERATURE = float(data['OUTER_CLOUD_SHIELD_TEMPERATURE'])
            self.INNER_CLOUD_SHIELD_TEMPERATURE = float(data['INNER_CLOUD_SHIELD_TEMPERATURE'])
            if (self.INNER_CLOUD_SHIELD_TEMPERATURE > self.OUTER_CLOUD_SHIELD_TEMPERATURE):
        		print "Bad cloud shield temperatures. Check the config file. Now using default variables..."
        		return False

            self.MINIMUM_DURATION = float(data['MINIMUM_DURATION'])
            self.MAXIMUM_DURATION = float(data['MAXIMUM_DURATION'])
            if (self.MINIMUM_DURATION > self.MAXIMUM_DURATION):
                print "Bad MIN and MAX DURATION inputs! Check the config file. Now using default variables..."
                return False;            
            self.DIRS = data['DIRS']

            # Checks that inputs are ok
            try:
                if not os.path.exists(self.DIRS['TRMMdirName']):
                    print "Error: TRMM invalid path!"
                    self.DIRS['TRMMdirName'] = raw_input("> Please enter the location to the raw TRMM netCDF files: \n")
            except:
                pass
            
            try:
                if not os.path.exists(self.DIRS['CEoriDirName']):
                    print "Error! MERG invalid path!"
                    self.DIRS['CEoriDirName'] = raw_input("> Please enter the directory to the MERG netCDF files: \n")
            except:
                print "..."   

            self.startDateTime = data['startDateTime']
            self.endDateTime = data['endDateTime']
            #check validity of time
            while utils.valid_date(self.startDateTime) != True:
                print "Invalid time entered for startDateTime!"
    
            while utils.valid_date(self.endDateTime) != True:
                print "Invalid time entered for endDateTime!"
            if (float(self.startDateTime) > float(self.endDateTime)):
                print "Bad start and end times input! Check the config file. Now using default variables..."
                return False   
                
            #check if all the files exisits in the MERG and TRMM directories entered
            test,_ = iomethods.check_for_files(self.DIRS['TRMMdirName'], self.startDateTime, self.endDateTime, 3, 'hour')
            if test == False:
                print "Error with files in the TRMM directory entered. Please check your files before restarting. "
                return

            test,self.filelist = iomethods.check_for_files(self.DIRS['TRMMdirName'], self.startDateTime, self.endDateTime, 3, 'hour')

        except ValueError:
            print "Bad config file, please check that inputs are float values. Now using default variables..."

        return True
        
class GraphVariables(object):
    def __init__(self):
        self.edgeWeight = [1,2,3]                   # weights for the graph edges
        self.CLOUD_ELEMENT_GRAPH = nx.DiGraph()     # graph obj for the CEs meeting criteria
        self.PRUNED_GRAPH = nx.DiGraph()            # graph meeting the CC criteria
        
def define_user_variables(useJSON=False):
    '''
    Purpose:: to actually assign the data according to a JSON file

    Inputs:: useJSON: a boolean indicating if the user has values in a JSON file

    Outputs:: a user variables instance with data initialized according to the userfile or defaults

    Assumptions::
    '''

    userVars = UserVariables()
    if useJSON: 
        if userVars.setupJSON(): 
            return userVars
        else:
            userVars2 = UserVariables() 
            return userVars2

    # we didn't want to use the config file, just return
    return userVars

def define_graph_variables():
    '''
    Purpose:: to assign the data for the GraphVariables

    Inputs:: None

    Outputs:: a Graphvariables instance 

    Assumptions::
    '''
    graphVars = GraphVariables()
    return graphVars


