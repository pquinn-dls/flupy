# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 19:15:05 2018

@author: pq67
"""

import numpy as np
from scipy import ndimage
from flupy.io.exchange_mantis import write_exchange_dataset_from_stack
import h5py

def rigid_shift_image(param,moving_img):
    _img = moving_img.copy()
    shape = _img.shape
    nx = shape[0]
    ny = shape[1]
    outofboundariesval = np.sum(_img)/float(nx*ny)             
    _img = ndimage.shift(_img,param,mode='constant', cval=outofboundariesval)
    return _img

#Apply image registration
def apply_image_registration(image, xshift, yshift):
        
    shape = image.shape
    nx = shape[0]
    ny = shape[1]
        
    outofboundariesval = np.sum(image)/float(nx*ny)        
    shifted_img = ndimage.interpolation.shift(image,[xshift,yshift],
                                                   mode='constant', 
                                                  cval=outofboundariesval)
        
    return shifted_img
    

def apply_alignment_to_stack(stack,xshifts,yshifts):
    
    nimages = stack.shape[2]
    for i in range(nimages):
        _img = stack[:,:,i]
        shape = _img.shape
        nx = shape[0]
        ny = shape[1]
        outofboundariesval = np.sum(_img)/float(nx*ny)             
        #_img = ndimage.rotate(_img,param[2],reshape=0)
        _img = ndimage.shift(_img,(xshifts[i],yshifts[i]),mode='constant', cval=outofboundariesval)
        stack[:,:,i] = _img
    return stack


def apply_alignment_to_file(infilename,outfilename, xshifts,yshifts):
    print 'Aligning the stack'
    fin = h5py.File(infilename,"r")
    data_stack = fin["exchange"]["data"][...]
    energies    = fin["exchange"]["energy"][...]
    fin.close()
    new_stack = apply_alignment_to_stack(data_stack,xshifts,yshifts)
    new_stack = crop_stack(new_stack,yshifts,xshifts)
    write_exchange_dataset_from_stack(outfilename,new_stack,energies)
    return new_stack,energies




def __convert_to_uint8bit(data, min_clip=None, max_clip=None, offset=None):
    '''
    Function converts input data to np.uint8 safely. Similar to the deprecated
    scipy.misc.bytescale.
    :param data. The ND input data as numpy array.
    :param min_clip. The minimum clip value. If None then data.min() is used.
                                                                Default:None.
    :param min_clip. The maximum clip value. If None then data.max() is used.
                                                                Default:None.
    :param offset. The offset from 0. If None then data.min() is used.
                                                                Default:None
    '''
    min_clip = data.min() if min_clip is None else min_clip
    max_clip = data.max() if max_clip is None else max_clip
    offset = min_clip if offset is None else offset
    intensity_range = max_clip - min_clip
    output_levels = 2 ** 8
    return ((data - offset) * (output_levels/intensity_range)).astype(np.uint8)


def crop_registed_images(images, min_rowshift, max_rowshift, min_colshift, max_colshift):
                
    n_rows  = images.shape[0]
    n_cols  = images.shape[1]
    # if the image is moved to the right (positive) we need to crop the left side 
    colleft = int(-1*(np.ceil(min_colshift)))
    if colleft < 0:
        colleft = 0
    # if the image is moved to the left (negative) we need to crop the right side 
    colright = int(np.floor(n_cols-max_colshift))
    if colright>(n_cols):
        colright = n_cols
        
    rowbottom = int(-1*(np.ceil(min_rowshift)))
    if rowbottom <0:
        rowbottom = 0
    rowtop = int(np.floor(n_rows-max_rowshift))
    if rowtop > (n_rows):
        rowtop = n_rows
    #print left,xright,ybottom,ytop
    print "Cropped range",rowbottom,rowtop,colleft,colright, min_rowshift, max_rowshift, min_colshift, max_colshift
    cropped_stack = images[ rowbottom:rowtop,colleft:colright, :]
    
    return cropped_stack #, xleft, xright, ybottom, ytop

def crop_stack(images, yshifts,xshifts):
    
    n_y, n_x = images.shape[0],images.shape[1]
    x_min = np.amin(xshifts[:])
    x_max = np.amax(xshifts[:])
    y_min = np.amin(yshifts[:])
    y_max = np.amax(yshifts[:])
    x1, y1 = int(np.round(y_min)), int(np.round(x_min))
    x2, y2 = int(np.round(y_max+0.5)), int(np.round(x_max+0.5))
    
    top1, bot1 = max(0, y1), min(n_y, n_y+y1)
    lef1, rig1 = max(0, x1), min(n_x, n_x+x1)
    top2, bot2 = max(0, y2), min(n_y, n_y+y2)
    lef2, rig2 = max(0, x2), min(n_x, n_x+x2)
    top = max(top1,top2)
    bot = min(bot1,bot2)
    lef = max(lef1,lef2)
    rig = min(rig1,rig2)
    print "crop1",top,bot,lef,rig,top1,top2,y1,y2
    print "crop2",y_min,y_max,x_min,x_max,n_y,n_x
    
    return images[top:bot,lef:rig,:]
