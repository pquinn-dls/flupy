import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import MaxNLocator


def scatter_histogram_plot(x,y,xlabel,ylabel,bins=50,filename=None):
    nullfmt   = NullFormatter()         # no labels
    
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    
    # start with a rectangular Figure
    plt.figure(1, figsize=(8,8))
    
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    # the scatter plot:
    pearR = np.corrcoef(x,y)[1,0]
    # least squares from:
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html
    A = np.vstack([x,np.ones(len(x))]).T
    m,c = np.linalg.lstsq(A,np.array(y))[0]

    axScatter.scatter(x, y)
    #axScatter.ylabel
    #axScatter.ylabel(ylabel)
    axScatter.plot(x,x*m+c,label="Fit %6.2e x + %6.2e\nPearson r = %6.2e"%(m,c,pearR))
    
    
    # now determine nice limits by hand:
    
#    xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
    xmax = np.max(x)
    ymax = np.max(y)
    xmin = np.min(x)
    ymin = np.min(y)
    
    xbinwidth = (xmax-xmin)/bins
    ybinwidth = (ymax-ymin)/bins
    
    
    axScatter.set_xlim( (xmin, xmax) )
    axScatter.set_ylim( (ymin, ymax) )
    
    xbins = np.arange(xmin+ xbinwidth, xmax + xbinwidth, xbinwidth)
    ybins = np.arange(ymin+ ybinwidth, ymax + ybinwidth, ybinwidth)

 #   print xbinwidth,xmin,xmax,xbins
 #   print ybinwidth,ymin,ymax,ybins
    axHistx.hist(x, bins=xbins)
    axHisty.hist(y, bins=ybins, orientation='horizontal')
    
    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )
#    axHistx.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    axHisty.xaxis.set_major_locator(MaxNLocator(prune='lower'))
#    axHisty.set_xscale("log")
#    axHistx.set_xscale("log")
    plt.setp(axHisty.get_xticklabels(), rotation=270, horizontalalignment='right')
    axScatter.legend(loc=0)
    axScatter.set_xlabel(xlabel)
    axScatter.set_ylabel(ylabel)
    if(filename !=None):
        plt.savefig(filename,bbox_inches='tight',dpi=100)
    else:
        plt.show()
    
    
#x = np.random.randn(1000)
#y = np.random.randn(1000)*500
#scatter_histogram_plot(x,y,"cr","co")

