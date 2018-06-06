from os import path
from flupy.utils.ConfigDict import ConfigDict

# This gets hooked to a general XRFDataset class
    
def loadDefaultConfigFiles(self):
    """
    Read in a configuration file using configdict
    Reads the dicts from the FXRFdata directory in pythonpath

    """
    dictlist=['FIT_PREFERENCES.DICT','DETECTOR.DICT','MATERIALS.DICT','ELEMENTS.DICT','OUTPUT_PREFERENCES.DICT','EXPT.DICT']
    paramdict = ConfigDict()
    this_dir, _this_filename = path.split(__file__)
    dictpath = path.join(this_dir, "FXRFdata")
    for dictname in dictlist:
        dictfile= path.join(dictpath,dictname)
        paramdict.read(dictfile)
    return paramdict

# Hook all these methods to XRFDataset

