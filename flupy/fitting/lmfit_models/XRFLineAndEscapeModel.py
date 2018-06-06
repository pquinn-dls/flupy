    
class XRFLineAndEscapeModel(Model):
    
    def __init__(self, paramdict,func,*args,**kwargs):
        """

        An XRF compton/scatter model
        LMFIT offers a range of options for coupling parameters and models but in order to be able to do 
        this it needs to manipulate parameters - in particular string expressions and this leads 
        big slowdowns if you've a large number of models or coupled parameters.		
        To speed up the calculation rather than create the XRF spectrum from a sum of Gaussian Models
        we create a single model with many parameters.
        The basic lmfit model works as follows:
        parse the arguments of the input function.
        Build a list of param names from these arguments
        Apply paramhints to these param_names or make parameters from the param names
  
        This class will pass a basic line

        """
        # Pass an empty function to model
        super(XRFScatterModel, self).__init__(None, *args, **kwargs)
        # Arguments need to be converted to  
        # a list of elements (the peaks -> translating to peak areas) + the line type arguments 
        # The basic model just assigns  
        # self.func = function
        # parses the function parameters to
        #
        self.paramdict    = paramdict
        self.extract_parameters()        
        
        
    def extract_parameters(self):
        # pull some commonly used values from the metadata
        self.incident_energy        = self.paramdict["Acquisition_instrument"]["beam_energy"]
        self.detector_type          = self.paramdict["Acquisition_instrument"]["XRF"]["Detector"]["detector_type"]
        self.lowE                   = self.paramdict['XRFFit']["FitParams"]["fitted_energy_range_keV"][0]
        self.highE                  = self.paramdict['XRFFit']["FitParams"]["fitted_energy_range_keV"][1]
        self.include_escape         = self.paramdict['XRFFit']["FitParams"]["include_escape"]
        

    def eval(self, params=None, **kwargs):
        """Evaluate the model with supplied parameters and keyword arguments.
        Parameters
        -----------
        params : Parameters, optional
            Parameters to use in Model.
        **kwargs : optional
            Additional keyword arguments to pass to model function.
        Returns
        -------
        numpy.ndarray
            Value of model given the parameters and other arguments.
        Notes
        -----
        1. if `params` is None, the values for all parameters are
        expected to be provided as keyword arguments.  If `params` is
        given, and a keyword argument for a parameter value is also given,
        the keyword argument will be used.
        2. all non-parameter arguments for the model function, **including
        all the independent variables** will need to be passed in using
        keyword arguments.
        """
        
        total=np.zeros_like(x)
        args = self.make_funcargs(params, kwargs)
        total = total + self.func(**args)
        peakcenternames = ["beam_energy","energy","centre","center"]
        for pname in args.keys():
        pname = model._strip_prefix(pname)
        if pname in peakcenternames:

        # now generate the escape peak:
        if self.include_escape:
            args["center"] = args["center"] - self.escape_energy
            args["area"]   = args["area"]*self.escape_factor
            total = total + self.func(**args)
        
        return total
    
