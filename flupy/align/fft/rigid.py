# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 20:42:44 2018

@author: pq67
"""
import h5py
import numpy as np
from skimage.feature import register_translation
import cv2
from flupy.align.utils import *

def multi_pass_alignment(filename, no_of_passes,reference=-1,nsmooth=2,imagelist=None):
    print 'Aligning the stack'
    fin = h5py.File(filename,"r")

    reference = reference
    data_stack = fin["exchange"]["data"][...]
    energies    = fin["exchange"]["energy"][...]
    if(imagelist):
        data_stack = data_stack[:,:,imagelist]
        energies    = energies[imagelist]
        
    print "Image shape",data_stack.shape
    reference_image = data_stack[:,:,reference]
    reference_image = cv2.blur(reference_image,(nsmooth,nsmooth))
    no_images = data_stack.shape[2]
    
    xshifts = np.zeros((no_images))
    yshifts = np.zeros((no_images))
    dxshifts = np.zeros((no_of_passes,no_images))
    dyshifts = np.zeros((no_of_passes,no_images))
    
    for j in range(no_of_passes):
        for i in range(no_images):
            img2 = data_stack[:,:,i] 
            img2 = cv2.blur(img2,(nsmooth,nsmooth))
            shifts,temp, ccorr = register_translation(reference_image,img2,upsample_factor=10,space="real")
            dxshifts[j,i] = shifts[0]
            dyshifts[j,i] = shifts[1]

        #Apply shifts
        xshifts[:] = xshifts[:]+dxshifts[j,:]
        yshifts[:] = yshifts[:]+dyshifts[j,:]
        
        for i in range(no_images):
            img = data_stack[:,:,i]
            if (abs(xshifts[i])>0.02) or (abs(yshifts[i])>0.02):
                shifted_img = apply_image_registration(img, dxshifts[j,i], dyshifts[j,i])
                data_stack[:,:,i] = shifted_img
        
        row_min = np.amin(dyshifts[j,:])
        row_max = np.amax(dyshifts[j,:])
        col_min = np.amin(dxshifts[j,:])
        col_max = np.amax(dxshifts[j,:])
        #data_stack =crop_registed_images(data_stack,col_min,col_max,row_min,row_max)
        #data_stack =crop_stack(data_stack,yshifts,xshifts)
        reference_image = data_stack[:,:,reference]
        reference_image = cv2.blur(reference_image,(nsmooth,nsmooth))
    fin.close()
    return data_stack,energies,xshifts,yshifts,dxshifts,dyshifts
