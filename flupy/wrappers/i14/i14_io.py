# -*- coding: utf-8 -*-
"""
Created on Thu Jun 07 16:40:20 2018

@author: pq67
"""

import h5py
import logging
import dask.array as da



def _hspy_checker(filename):
    if h5py.is_hdf5(filename):
        hdf5_file = h5py.File(filename, mode='r')
        for attr in hdf5_file.attrs:
            if 'file_format' in attr:
                if hdf5_file.attrs['file_format'] == 'HyperSpy':
                    return(True)
    return(False)


def _load_lazy_i14_file(filename, chunk_size=(16, 16)):
    f = h5py.File(filename, mode='r')
    if 'entry1/xspress3detector' in f:
        data = f['/fpd_expt/fpd_data/data']
        if len(data.shape) == 4:
            chunks = (
                    chunk_size[0], chunk_size[1],
                    1, data.shape[-2], data.shape[-1])
            data_lazy = da.from_array(data, chunks=chunks)[:, :, 0, :, :]
            s_out = LazyPixelatedSTEM(data_lazy)
        elif len(data.shape) == 3:
            s_out = _load_fpd_emd_file(filename, lazy=True)
        else:
            raise IOError(
                "I14 XRF dataset does not have correct dimensions")
        return(s_out)
    else:
        raise IOError("I14 XRF dataset not found")


def _load_other_file(filename, lazy=False):
    s = load(filename, lazy=lazy)
    s_new = signal_to_pixelated_stem(s)
    return s_new


def signal_to_pixelated_stem(s):
    """Make a PixelatedSTEM object from a HyperSpy signal.

    This will retain both the axes information and the metadata.
    If the signal is lazy, the function will return LazyPixelatedSTEM.

    Parameters
    ----------
    s : HyperSpy signal
        Should work for any HyperSpy signal.

    Returns
    -------
    pixelated_stem_signal : PixelatedSTEM or LazyPixelatedSTEM object

    Examples
    --------
    >>> import numpy as np
    >>> import hyperspy.api as hs
    >>> import fpd_data_processing.api as fp
    >>> s = hs.signals.Signal2D(np.random.random((8, 11, 21, 13)))
    >>> s.metadata.General.title = "test dataset"
    >>> s
    <Signal2D, title: test dataset, dimensions: (11, 8|13, 21)>
    >>> from fpd_data_processing.io_tools import signal_to_pixelated_stem
    >>> s_new = signal_to_pixelated_stem(s)
    >>> s_new
    <PixelatedSTEM, title: test dataset, dimensions: (11, 8|13, 21)>

    """
    # Sorting axes as a function of its index
    axes_list = [x for _, x in sorted(s.axes_manager.as_dictionary().items())]
    metadata = s.metadata.as_dictionary()
    if s._lazy:
        s_new = LazyPixelatedSTEM(s.data, axes=axes_list, metadata=metadata)
    else:
        s_new = PixelatedSTEM(s.data, axes=axes_list, metadata=metadata)
    return s_new


def load_xrf_signal(
        filename, lazy=False, chunk_size=(16, 16),
        navigation_signal=None):
    """
    Parameters
    ----------
    filename : string
    lazy : bool, default False
    chunk_size : tuple, default (16, 16)
        Used if Lazy is True. Sets the chunk size of the signal in the
        navigation dimension. Higher number will potentially make the
        calculations be faster, but use more memory.
    navigation_signal : Signal1D

    """
    
    if _fpd_checker(filename, attr_substring='fpd_version'):
        if lazy:
            s = _load_lazy_fpd_file(filename, chunk_size=chunk_size)
        else:
            s = _load_fpd_emd_file(filename)
    elif _hspy_checker(filename, attr_substring='HyperSpy'):
        s = _load_other_file(filename, lazy=lazy)
    else:
        # Attempt to load non-fpd and non-HyperSpy signal
        s = _load_other_file(filename, lazy=lazy)
    if navigation_signal is None:
        try:
            s_nav = _load_fpd_sum_im(filename)
            s.navigation_signal = s_nav
        except IOError:
            logging.debug("Nav signal not found in {0}".format(filename))
            s.navigation_signal = None
        except ValueError:
            logging.debug("Nav signal in {0}: wrong shape".format(filename))
            s.navigation_signal = None
    else:
        nav_im_shape = navigation_signal.axes_manager.signal_shape
        nav_ax_shape = s.axes_manager.navigation_shape
        if nav_im_shape == nav_ax_shape:
            s.navigation_signal = navigation_signal
        else:
            raise ValueError(
                    "navigation_signal does not have the same shape ({0}) as "
                    "the signal's navigation shape ({1})".format(
                        nav_im_shape, nav_ax_shape))
    return s

