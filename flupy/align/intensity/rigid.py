# -*- coding: utf-8 -*-
from scipy import ndimage
import numpy as np
from scipy import optimize
from .metrics import mutual_information,mutual_information_2d
from flupy.flupy.align.utils import rigid_shift_image

def MIErrFunc(param,ref_img,moving_img,bins=20):
    
    # Perform a translational shift
    _img =rigid_shift_image(param,moving_img)
    # measure the mutual information between the reference and moving image
    return mutual_information(ref_img,_img,bins)


def MIErrFunc2(param,ref_img,moving_img,bins=20):
    
    # Perform a translational shift
    _img =rigid_shift_image(param,moving_img)
    # measure the mutual information between the reference and moving image
    return mutual_information_2d(ref_img,_img,bins=bins)



def rigidregistration(img,ximg,errFunc,runCoarse=True,opt='global'):
    # Rough search...
    # rough estimate of number of bins based on the starting image size...
    xsize,ysize = img.shape
    bins = int(np.sqrt(xsize*ysize/5.))

    xsize = int(0.4*xsize)
    ysize = int(0.4*ysize)

    if opt=='global':
        bounds=[(-xsize,xsize),(-ysize,ysize)]
        param = optimize.differential_evolution(errFunc,bounds,args=(img,ximg,bins))
        param = param.x
        _img = rigid_shift_image(param,ximg)
    else:
        x_range = np.array(xrange(-xsize,xsize))
        y_range = np.array(xrange(-ysize,ysize))
        bestx = 0
        besty = 0
        besterr = 1e10
        for xpos in x_range:
            for ypos in y_range:
                err = errFunc([xpos,ypos],img,ximg)
                if err < besterr:
                    besterr=err
                    bestx = xpos
                    besty = ypos
        param = [bestx,besty]
        param = optimize.fmin(errFunc,param,args=(img,ximg,bins))
        
    return (_img,param)