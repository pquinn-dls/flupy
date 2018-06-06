import os
import numpy as np
import re
import h5py
import time
from flupy.algorithms.xrf_calculations.solid_angle import solid_angle_fraction

"""

Load the data...
In this case from a DLS NEXUS file...

Returns a dictionary object with

"data" - xrf data  (rows x cols x mca_length)
"average"  - average xrf spectrum calculated from the data
"rows" - number of rows
"cols" - number of columns    
"i0"  beam intensity for normalization if available
"xy"  xy grid positions if available 

+Metadata....

Basically metadata is stored in the form of the configuration data
This is so that it can be added to or overwrite the configuration
depending on the data...
For example, the incident energy, the detector used, the detector distance,
calibration data, collection times, may all be stored in the hdf or nexus file
but used to help the fitting process.  


"""

def xrfresultwrite(infile,paramdict,fitdict,datadict,matrixdict, overwrite=True,writergb=True,fullresults=False):
    
    #
    # write out the rgb files...
    #    
    # create_chisq, elements,back,scatter rgb file
    #
    # Store the same data in a nexus file...
    #
    # Create file to output data. This allows for larger files to be processed
    # by the program.
    #
    print "Writing nexus and rgb file"
    outfile = os.path.splitext(os.path.basename(infile))[0]
    dirname = os.path.dirname(infile)
    dirname = os.path.join(dirname,"IMAGES")
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    outfile=outfile + '_'
    if("output_dir" in paramdict["Output"]):
        if(not paramdict["Output"]["output_dir"]):
            output_dir=""
        else:    
            output_dir=paramdict["Output"]["output_dir"]
    else:
        output_dir=""
    
    timestr = time.strftime("%Y%m%d_%H%M%S",time.localtime())
    uoutfile = os.path.join(output_dir,outfile + timestr + '.hdf5')
    print "nexus output file",uoutfile
    dataOut = h5py.File(uoutfile,"w")

    inpgrp = dataOut.create_group("entry1")
    inpgrp.attrs["XRFfile"]=np.string_(infile)
    fitmat = inpgrp.create_group("InputMatrices")
    
    
    DB= matrixdict["Total_Matrix"]
    Ap= matrixdict["Weighted_Total_Matrix"]  
    DB_constraints=matrixdict["Total_Constraints_Matrix"] 
    rweight= matrixdict["Weighting"] 
    
    DB_scatter= matrixdict["Scatter_Matrix"]  
    DB_scatter_constraint=matrixdict["Scatter_Constraints"]  
    DB_scatter_descript=matrixdict["Scatter_Description"]  
    DB_scatter_areas=matrixdict["Scatter_Areas"]  
    nscatt_params=matrixdict["Scatter_Num_curves"]  
       
    DB_back=matrixdict["Background_Matrix"]       
    DB_back_constraint=matrixdict["Background_Constraints"]  
    DB_back_areas=matrixdict["Background_Areas"]  
    DB_back_descript=matrixdict["Background_Description"]  
    nback_params=matrixdict["Background_Num_curves"]  
    
    DB_el=matrixdict["XRF_Matrix"]      
    DB_el_constraint=matrixdict["XRF_Constraints"]  
    DB_el_areas=matrixdict["XRF_Areas"]  
    DB_el_cs=matrixdict["XRF_CrossSec"]      
    DB_el_descript=matrixdict["XRF_Description"]  
    nel_params=matrixdict["XRF_Num_curves"]  
    if(DB_scatter!=None):
        matrixdict["Descriptions"] = DB_el_descript + DB_back_descript + DB_scatter_descript
    else:
        matrixdict["Descriptions"] = DB_el_descript + DB_back_descript    
        
        
    fitmat.create_dataset("RawMatrix",data=DB)
    fitmat.create_dataset("ElementAreas",data=DB_el_areas)            
    fitmat.create_dataset("BackgroundAreas",data=DB_back_areas)                
    if(nscatt_params>0):
        fitmat.create_dataset("ScatterAreas",data=DB_scatter_areas)  
    mylist=[]
    for el in DB_el_descript:    
        mylist.append(el[2]+"_"+el[1])
    for i in range(nback_params):
        mylist.append("background"+str(i))
    if(nscatt_params>0):    
        for i,scatt in enumerate(DB_scatter_descript):
            mylist.append(scatt[0]+str(i)+"_"+scatt[1])
    
    stringlist=np.array(mylist)
    fitmat.create_dataset("ParamLists",data=stringlist)  
    
    outgrp = inpgrp.create_group("OutputResults")
    fitres = outgrp.create_group("Maps")
    

    pdq1=np.multiply(fitdict["parameters"][:,:,:nel_params],np.sum(DB_el_areas,axis=1)/np.sum(DB_el_areas>0,axis=1))
    correction = paramdict["Experiment"]["photon_flux"]
    print "photon flux",correction
    if("Matrix" in paramdict):
        density=paramdict["Matrix"]["Density"]
        thickness=paramdict["Matrix"]["Thickness"]
        correction = correction*density*thickness
        print "Matrix correction"
    prefactor = paramdict["Output"]["unit_prefactor"] 
    unit_label=paramdict["Output"]["unit_label"]
    solid_angle = solid_angle_fraction(paramdict)
    #solid_angle=1.0
    collectionTime=datadict["Experiment"]["collection_time"] 
    rows =datadict["rows"] 
    cols =datadict["cols"]     
    correction = correction*prefactor*collectionTime*solid_angle
    print "correction!",correction,solid_angle,collectionTime
    pdq1=pdq1/correction
    # clean up rounding errors which produce really small values....
    
    pdq1= np.where(pdq1<0.0,1.0e-9,pdq1)
    windowlist=[]
    for i,el in enumerate(DB_el_descript):
        fitres.create_dataset(el[2]+"_"+el[1]+unit_label,data=pdq1[:,:,i])
        windowlist.append(el[2]+"_"+el[1]+unit_label)

    mest = fitdict["parameters"]
    chi_sq=fitdict["chisq"]
    if(nscatt_params>0):
        pdq2=np.tensordot(DB[:,nel_params+nback_params:nel_params+nback_params+nscatt_params], mest[:,:,nel_params+nback_params:nel_params+nback_params+nscatt_params], (1,2))
        pdq2=np.sum(pdq2,axis=0)
        fitres.create_dataset("scatter",data=pdq2)
        windowlist.append("scatter")
    

    pdq3=np.tensordot(DB[:,nel_params:nel_params+nback_params], mest[:,:,nel_params:nel_params+nback_params], (1,2))
    pdq3=np.sum(pdq3,axis=0)

    fitres.create_dataset("background",data=pdq3)
    windowlist.append("background")


    fitres.create_dataset('chisq',data=chi_sq)  
    windowlist.append("chisq")
    if(nscatt_params>0):
        pdq1 = np.dstack((pdq1,pdq2,pdq3,chi_sq))
    else:
        pdq1 = np.dstack((pdq1,pdq3,chi_sq))
    fitres = outgrp.create_group("FitParameters")
    fitres.create_dataset('parameters',data=mest)  
    dataOut.close()
    print "Finished writing nexus file"

    if(writergb):
        rgbfile = os.path.join(output_dir,outfile + ".rgb")
        print "rgb output file",rgbfile
        fout=open(rgbfile,"w")
        mystr='row  column  '
        for k in range(len(windowlist)):
            mystr=mystr+windowlist[k]
            if(k<(len(windowlist)-1) ):
                mystr=mystr+"  "
        print >> fout, mystr
        for row in range(rows):
            for col in range(cols):
                stp=" ".join(map(str, pdq1[row,col].astype('str')))
                print >>fout,row,col,stp    
        # close the files...
        fout.close()            
        print "Finished writing rgb file :",rgbfile

    return rgbfile,uoutfile
