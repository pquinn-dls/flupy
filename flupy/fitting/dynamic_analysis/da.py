#
#
# Dynamic Analysis methods
#
#
# Metadata format
#
#├── DynamicAnalysis
#│   ├── Weighted  * Boolean
#│   ├── energy    * array of points
#│   ├── a_inv     * using DA paper naming convention
#│   ├── b
#│   ├── bb
#│   ├── a_inv*b
#│   ├── lines     * list of line names (elements...,pileup,scatter,compton)
#├── Sample
#│   ├── description  
#│   ├── elements
#│   ├── xray_lines
#│   ├── composition
#│   ├── concentrations        
#│   ├── matrix        
#│   └── multilayer
#
#
#
#
def loadDAMatrix(filename):
    """
	Loads a Dynamic Analysis Matrix from a file.
	See:
	
	The metadata structure is: 
		DA
	├── Matrix	
	│   ├── Weighted  :Boolean
	│   ├── energy    :Numpy array 
	│   ├── a_inv     :Numpy array - arrays namesd using DA paper naming convention
	│   ├── b         :Numpy array 
	│   ├── bb        :Numpy array
	│   ├── a_inv*b   :Numpy array 
	│   ├── lines     :String list - corresponding to the DA spectra or rows in the array  (e.g. "Fe", "pileup","scatter","compton")
	├── Yield   *if present 	
	│   ├── energy                *  array of points
	│   ├── tranmssion_filters    *  array of points evaluated at energy`  
	│   ├── detector_efficiency   *  array of points evaluated for each line/element 
	│   ├── photon_flux            
	│   ├── calibration_per_line   
	├── Sample
	│   ├── description    :string  
	│   ├── elements       :string list  e.g  ["Fe","Co"]
	│   ├── xray_lines     :string list  e.g. ["Fe_K"] or more specific descriptions ["Fe_Ka","Fe_Kb"]
	│   ├── composition    :
	│   ├── concentrations :        
	│   ├── matrix         : 
	│   └── multilayer     : 
    Parameters
    ----------
    filename : string
        HDF5 file from which to load the DA matrix

    Returns
    -------
    DAdict : dict
        A python dictionary containing the data above
  
		
    """
	

def saveDAMatrix(filename,DAdict):


def saveDA(filename,DAdict):


def loadYieldFile(filename,DAdict):


def saveYieldFile(filename,DAdict):


def loadParameterDict(dictFileLocation=None):


def createDAMatrixFromModel(model,createYield=True):


def createDAMatrixFromlmfit(params,model):
    """
    
    """
    
    #
    # rather than trying to overwrite the lmfit model we will basically cycle 
    # through
    #
    



def createYieldMatrixFromParameters(paramdict):


def solveSpectrum(DAdict,spectrum):
    """
    
    
    
    
    
    """
    return 
