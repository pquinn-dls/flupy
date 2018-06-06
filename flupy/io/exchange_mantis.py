# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 19:12:45 2018

@author: pq67
"""
import numpy as np
import h5py

def write_exchange_dataset_from_stack(file_name,image_stack,energies):
    #And adds the necessary groups, datasets and attributes(Scientific Data Exchange standard HDF5 format data)
    xnumber = np.arange(image_stack.shape[0]*1.0)
    ynumber = np.arange(image_stack.shape[1]*1.0)
    enumber = np.arange(image_stack.shape[2]*1.0)
    #print xnumber.shape,ynumber.shape,xnumber,ynumber
    inumber = np.ones(image_stack.shape[2])
    comment=''
    f1 = h5py.File(file_name,'w')
    dset = f1.create_group("exchange")
    dset2 = dset.create_dataset("data",data=image_stack)
    dset2.attrs['axes'] = 'x:y'
    dset2.attrs['signal'] = 1
    dset3 = dset.create_dataset("energy",data=energies)
    dset3.attrs['units'] = 'eV'
    dset4 = dset.create_dataset("x",data=xnumber)
    dset5 = dset.create_dataset("y",data=ynumber)
    str_type = h5py.new_vlen(str)
    eset = f1.create_dataset("implements", shape=(1,), dtype=str_type)
    eset[:] = 'information:exchange:spectromicroscopy'
    fset = f1.create_group("information")
    fset2 = fset.create_dataset("comment", shape=(1,), dtype=str_type)
    fset2[:] = comment
    fset3 = fset.create_dataset("file_creation_datetime", shape=(1,), dtype=str_type)
    fset3[:] = "2012-07-11T09:15"
    fset3.attrs['file_creation_datetime'] = 'time'
    gset = f1.create_group("spectromicroscopy")
    gset2 = gset.create_group("normalization")
    gset3 = gset2.create_dataset("white_spectrum",data=inumber)
    gset4 = gset2.create_dataset("white_spectrum_energy", data=enumber)
    gset4.attrs['units'] = 'eV'
    hset = f1.create_dataset("version", shape=(1,), dtype=str_type)
    hset[:] = '1.0'
    f1.close()
   
