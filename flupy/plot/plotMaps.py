import numpy as np
#import pylab
#import matplotlib.pyplot as plt
#from matplotlib.transforms import BlendedGenericTransform
import pylab

def plotMaps(paramdict,fitdict,datadict,matrixdict,xp,yp):
        """
        
        A worker routine
        Plots the best and worst fit results at the end of the map
        Useful to check the worst one and see if the background or other paramaters
        need changing... 
        
        """
        titles = "Fit at Position",""
        mylist  = []

        labels=[]
        nel_params = matrixdict["XRF_Num_curves"]
        nback_params=matrixdict["Background_Num_curves"]
        nscatt_params=matrixdict["Scatter_Num_curves"]
        DB  = matrixdict["Total_Matrix"]
        descript_el=matrixdict["XRF_Description"]
        result_est=fitdict["parameters"]
        x=paramdict["FitParams"]["mca_energies_used"]
        irange=paramdict["FitParams"]["mca_channels_used"]
        data = datadict["data"]
        fig=pylab.figure()
        nplots = 100*len(mylist)+10
        myplots = []
        for i,fit in enumerate(mylist):
            labels=[]
            ax = pylab.subplot(nplots+i)
            back    = np.dot(DB[:,nel_params:nel_params+nback_params],result_est[fit[0],fit[1]][nel_params:nel_params+nback_params])
            result = np.dot(DB,result_est[fit[0],fit[1]])
            for j,el in enumerate(descript_el):
                if(el[1]=='xrf'):
                    a, = ax.plot(x, back+DB[:,j]*result_est[fit[0],fit[1]][j],label=el[2])
                    labels.append(el[2])
                    myplots.append(a)
                elif(el[1]=='pileup'):
                    a,=ax.plot(x,back+DB[:,j]*result_est[fit[0],fit[1]][j],label=el[2]+"_pileup")
                    labels.append(el[2]+"_pileup")
                    myplots.append(a)
                elif(el[1]=='escape'):
                    a,=ax.plot(x, back+DB[:,j]*result_est[fit[0],fit[1]][j],label=el[2]+"_escape")
                    labels.append(el[2]+"_escape")
                    myplots.append(a)        
            ax.set_title(titles[i])
            a,= ax.plot(x,result,linestyle='--',label="Fit")
            labels.append("Fit")
            
            myplots.append(a)        
            a,=ax.plot(x,back,label="background")
            labels.append("background")
            myplots.append(a)    
            if(nscatt_params>0):
                scatt   = np.dot(DB[:,nel_params+nback_params:nel_params+nback_params+nscatt_params],\
                                 result_est[fit[0],fit[1]][nel_params+nback_params:nel_params+nback_params+nscatt_params])
                a,=ax.plot(x, scatt,label="scatter")
                labels.append("scatter")
                myplots.append(a)
            a,=ax.plot(x,data[fit[0],fit[1]][irange],label="raw data")
            labels.append("raw data")
            myplots.append(a)

        _leg = pylab.figlegend(myplots,labels,loc=10,bbox_to_anchor=(0.9, 0.5),prop={'size':8})
        fig.subplots_adjust(right=0.8)
        pylab.show()
