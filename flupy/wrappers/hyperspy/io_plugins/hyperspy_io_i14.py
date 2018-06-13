# -*- coding: utf-8 -*-
"""
Created on Thu Jun 07 16:40:20 2018

@author: pq67
"""

import h5py
import logging
import dask.array as da
from hyperspy.io import load
from hyperspy.signals import Signal1D


def _hspy_checker(filename):
    if h5py.is_hdf5(filename):
        hdf5_file = h5py.File(filename, mode='r')
        for attr in hdf5_file.attrs:
            if 'file_format' in attr:
                if hdf5_file.attrs['file_format'] == 'HyperSpy':
                    return(True)
    return(False)


dataloc_i14=''
mapping_i14={}

dataloc_i20=''
mapping_i20={}



dataloc_i18_xmap=''
mapping_i18_xmap={}



def _load_lazy_dls_xrf_file(filename, chunk_size=(16, 16),dataloc,mapping):
    #
    # first load the metadata
    #
    metadata = loaddict()
    #
    # Now load the data
    #
    f = h5py.File(filename, mode='r')
    if dataloc in f:
        data = f[dataloc]
        chunks = (
                chunk_size[0], chunk_size[1], data.shape[-2], data.shape[-1])
        data_lazy = da.from_array(data, chunks=chunks)
        # so we have the data
        #
        # now extract some metadata from the nexus file....
        #
        #  
        dims = data_lazy.size
        
        axes = [{
            'size': dims[0],
                'index_in_array': 0,
                'name': 'X',
                'scale': 1,
                'offset': 0,
                'units': 'mm',
                'navigate': True,
                },
               {
            'size': dims[1],
                'index_in_array': 0,
                'name': 'Y',
                'scale': 1,
                'offset': 0,
                'units': 'mm',
                'navigate': True,
                },
               {
            'size': dims[2],
                'index_in_array': 0,
                'name': 'detector',
                'scale': 1,
                'offset': 0,
                'navigate': True,
                },
               {
            'size': dims[3],
                'index_in_array': 0,
                'name': 'channels',
                'scale': 0.01,
                'offset': 0.,
                'units': 'keV',
                'navigate': False,
                }]

        dictionary = {'data': data,
                      'axes': axes,
                      'metadata': metadata,
                      'mapping': mapping,
                      'original_metadata': original_metadata}

    return [dictionary, ]
        
        s_out = LazyXRFSpectrum(data_lazy)
        return(s_out)
    else:
        raise IOError("I14 XRF dataset not found")



def file_reader(filename, *args, **kwds):
    with open(filename, 'rt') as f:
        # Strip leading, empty lines
        line = str(f.readline())
        while line.strip() == '' and not f.closed:
            line = str(f.readline())
        try:
            date, version = line.split('\t')
        except ValueError:
            _bad_file(filename)
        if version.strip() != 'Digiheater 3.1':
            _bad_file(filename)
        calib = str(f.readline()).split('\t')
        str(f.readline())       # delta_t
        header_line = str(f.readline())
        try:
            R0, a, b, c = [float(v.split('=')[1]) for v in calib]
            date0 = datetime.strptime(date, "%d/%m/'%y %H:%M")
            date = '%s' % date0.date()
            time = '%s' % date0.time()
        except ValueError:
            _bad_file(filename)
        original_metadata = dict(R0=R0, a=a, b=b, c=c, date=date0,
                                 version=version)

        if header_line.strip() != (
                'sample\ttime\tTset[C]\tTmeas[C]\tRheat[ohm]\tVheat[V]\t'
                'Iheat[mA]\tPheat [mW]\tc'):
            _bad_file(filename)
        try:
            rawdata = np.loadtxt(f, converters={1: _cnv_time}, usecols=(1, 3),
                                 unpack=True)
        except ValueError:
            _bad_file(filename)

    times = rawdata[0]
    # Add a day worth of seconds to any values after a detected rollover
    # Hopefully unlikely that there is more than one, but we can handle it
    for rollover in 1 + np.where(np.diff(times) < 0)[0]:
        times[rollover:] += 60 * 60 * 24
    # Raw data is not necessarily grid aligned. Interpolate onto grid.
    dt = np.diff(times).mean()
    temp = rawdata[1]
    interp = scipy.interpolate.interp1d(times, temp, copy=False,
                                        assume_sorted=True, bounds_error=False)
    interp_axis = times[0] + dt * np.array(range(len(times)))
    temp_interp = interp(interp_axis)

    metadata = {'General': {'original_filename': os.path.split(filename)[1],
                            'date': date,
                            'time': time},
                "Signal": {'signal_type': "",
                           'quantity': "Temperature (Celsius)"}, }

    axes = [{
        'size': len(temp_interp),
            'index_in_array': 0,
            'name': 'X',
            'scale': dt,
            'offset': times[0],
            'units': 'mm',
            'navigate': True,
            }]

    dictionary = {'data': temp_interp,
                  'axes': axes,
                  'metadata': metadata,
                  'original_metadata': {'DENS_header': original_metadata},
                  }

    return [dictionary, ]


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


def load_i14_xrf_signal(
        filename, dictFolder='auto',lazy=False, chunk_size=(16, 16),
        navigation_signal=None):
    """
    
    A signal needs a dict containing:
        
            data : numpy array
               The signal data. It can be an array of any dimensions.
            axes : dictionary (optional)
                Dictionary to define the axes (see the
                documentation of the AxesManager class for more details).
            attributes : dictionary (optional)
                A dictionary whose items are stored as attributes.
            metadata : dictionary (optional)
                A dictionary containing a set of parameters
                that will to stores in the `metadata` attribute.
                Some parameters might be mandatory in some cases.
            original_metadata : dictionary (optional)
                A dictionary containing a set of parameters
                that will to stores in the `original_metadata` attribute. It
                typically contains all the parameters that has been
                imported from the original data file.


    Parameters
    ----------
    filename   : string
    dictfolder : string - 'Auto' or a directory containing the metadata
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

def load_dls_xrf_signal(
        filename, dataloc,position_loc=None,sum_detectors=True,
        dictFolder='auto',lazy=False, chunk_size=(16, 16),
        ):

    
    
def _load_lazy_fpd_file(filename, chunk_size=(16, 16)):
    f = h5py.File(filename, mode='r')
    if 'fpd_expt' in f:
        data = f['/fpd_expt/fpd_data/data']
        if len(data.shape) == 5:
            chunks = (
                    chunk_size[0], chunk_size[1],
                    1, data.shape[-2], data.shape[-1])
            data_lazy = da.from_array(data, chunks=chunks)[:, :, 0, :, :]
            s_out = LazyPixelatedSTEM(data_lazy)
        elif len(data.shape) == 4:
            s_out = _load_fpd_emd_file(filename, lazy=True)
        else:
            raise IOError(
                "Pixelated dataset does not have correct dimensions")
        return(s_out)
    else:
        raise IOError("Pixelated dataset not found")


def _load_fpd_sum_im(filename):
    f = h5py.File(filename, mode='r')
    if 'fpd_expt' in f:
        if len(f['/fpd_expt/fpd_sum_im/data'].shape) == 3:
            data = f['/fpd_expt/fpd_sum_im/data'][:, :, 0]
        elif len(f['/fpd_expt/fpd_sum_im/data'].shape) == 2:
            data = f['/fpd_expt/fpd_sum_im/data'][:, :]
        else:
            ValueError("fpd_sum_im does not have the correct dimensions")
        s = Signal2D(data)
        f.close()
        return(s)
    else:
        raise IOError("Pixelated dataset not found")


def _load_fpd_sum_dif(filename):
    f = h5py.File(filename, mode='r')
    if 'fpd_expt' in f:
        if len(f['/fpd_expt/fpd_sum_dif/data'].shape) == 3:
            data = f['fpd_expt/fpd_sum_dif/data'][0, :, :]
        elif len(f['/fpd_expt/fpd_sum_dif/data'].shape) == 2:
            data = f['fpd_expt/fpd_sum_dif/data'][:, :]
        else:
            ValueError("fpd_sum_dif does not have the correct dimensions")
        s = Signal2D(data)
        f.close()
        return(s)
    else:
        raise IOError("Pixelated dataset not found")


def _load_fpd_emd_file(filename, lazy=False):
    logging.basicConfig(level=logging.ERROR)
    s_list = load_with_reader(filename, reader=emd, lazy=lazy)
    logging.basicConfig(level=logging.WARNING)
    temp_s = None
    longest_dims = 0
    for s in s_list:
        if len(s.data.shape) > longest_dims:
            longest_dims = len(s.data.shape)
    if longest_dims < 4:
        raise ValueError("Pixelated dataset not found")
    for s in s_list:
        if len(s.data.shape) == longest_dims:
            temp_s = s
            break
    if longest_dims == 4:
        s = temp_s.transpose(signal_axes=(0, 1))
    elif longest_dims == 5:
        s = temp_s.isig[:, :, 0, :, :].transpose(signal_axes=(0, 1))
    else:
        raise Exception(
                "Pixelated dataset not found")
    s._lazy = lazy
    s_new = signal_to_pixelated_stem(s)
    return(s_new)


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


def load_fpd_signal(
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
    navigation_signal : Signal2D
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

