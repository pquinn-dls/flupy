from ConfigDict import ConfigDict
import os
#
# Routines to load and save python dictionaries 
# 
#


def loadDictFromFiles(path_to_config, config_file_names):
    """
    
	 Read in a configuration file using configdict
    Reads the dicts from the FXRFdata directory in pythonpath
	
	 Parameters:
	
	 path_to_config    :  string
	 Directory containing the config files
	
	 config_file_names : string list  
	 List containing file names to be loaded
	
	 Returns:
	
	 paramdict:  A dictionary containing the contents of the various 
	 dict files which were written using the configdict or configparser
	 methods and formats
	
    """

    dict_list=[]
    for localfile in config_file_names:
        fullpath = os.path.join(path_to_config, localfile)
        if os.path.isfile(fullpath):     
            dict_list.append(fullpath)
    if dict_list:
        paramdict = ConfigDict()
        paramdict.read(dict_list)
        return paramdict
    else:
        return None



def loadDictFromDir(path_to_config,extension='.dict'):
    
    """
    Read in a configuration file using configdict
    Reads the dicts from the FXRFdata directory in pythonpath
	
	 Parameters:
	 path_to_config    :  string
	 Directory containing the config files
	
    config_file_names : string list  
    List containing file names to be loaded
	
    Returns:
	
    paramdict:  A dictionary containing the contents of the various 
    dict files which were written using the configdict or configparser
    methods and formats
	
    """
    myext = extension
    dict_list=[]
    for localfile in os.listdir(path_to_config):
        if localfile.endswith(myext):
            dict_list.append(os.path.join(path_to_config, localfile))
    print dict_list,len(dict_list)
    if dict_list:
        paramdict = ConfigDict()
        paramdict.read(dict_list)
        return paramdict
    else:
        return None


def loadDictFromDefault(extension='.dict'):
    
    """
    Read in a configuration file using configdict
    Reads the dicts from the FXRFdata directory in pythonpath
	
	 Parameters:
	 path_to_config    :  string
	 Directory containing the config files
	
    config_file_names : string list  
    List containing file names to be loaded
	
    Returns:
	
    paramdict:  A dictionary containing the contents of the various 
    dict files which were written using the configdict or configparser
    methods and formats
	
    """
    myext = extension
    dict_list=[]
    this_dir, _this_filename = os.path.split(__file__)
    this_dir = os.path.join(this_dir,"../db/data/")
    path_to_config =os.path.abspath(this_dir)

    for localfile in os.listdir(path_to_config):
        if localfile.endswith(myext):
            dict_list.append(os.path.join(path_to_config, localfile))
    
    if dict_list:
        paramdict = ConfigDict()
        paramdict.read(dict_list)
        return paramdict
    else:
        return None
	

def loadFileFromDefault(localfile):
    
    """
    Read in a configuration file using configdict
    Reads the dicts from the FXRFdata directory in pythonpath
	
	 Parameters:
	 path_to_config    :  string
	 Directory containing the config files
	
    config_file_names : string list  
    List containing file names to be loaded
	
    Returns:
	
    paramdict:  A dictionary containing the contents of the various 
    dict files which were written using the configdict or configparser
    methods and formats
	
    """
    dict_list=[]
    this_dir, _this_filename = os.path.split(__file__)
    this_dir = os.path.join(this_dir,"../db/data/")
    path_to_config =os.path.abspath(this_dir)
    if not localfile.endswith('.dict'):
        localfile + ".dict"
    dict_list.append(os.path.join(path_to_config, localfile))
    
    if dict_list:
        paramdict = ConfigDict()
        paramdict.read(dict_list)
        return paramdict
    else:
        return None



def saveDictByKeys(userdict,path_to_config):
    
    """
    
	 Save the dictionary to a set of files
    The dictionary is split using the top level keys i.e. dict.keys()
    A seperate .DICT file is created for each key.
 	
	 Parameters:
    userdict          :  python dictionary	
	 path_to_config    :  string path of folder to save to dict files to
	
    """
    outputdict = ConfigDict(userdict)
    outputdict.write_sections(path_to_config)

def saveDictToFile(userdict,filename):
    
    """
    Save a dictionary to a configuration file
	 Parameters:
    userdict        :  python dictionary
	 filename    :  string - full file path name
	
	
    """
    outputdict = ConfigDict(userdict)
    outputdict.write(filename)
	
	
	
	