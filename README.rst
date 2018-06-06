flupy
========

.. image:: https://travis-ci.org/lmfit/lmfit-py.png
   :target: https://travis-ci.org/lmfit/lmfit-py

.. image:: https://zenodo.org/badge/4185/lmfit/lmfit-py.svg
   :target: https://zenodo.org/badge/latestdoi/4185/lmfit/lmfit-py


Overview
---------

flupy provides various tools for x-ray fluorescence (XRF) fitting and analysis
The package is based around the needs of the I14 beamline to address:

    Basic x-ray calculations
	Non-linear XRF Fitting using LMfit - allow for custom and user peaks 
	Linear fitting using dynamic analysis and nnls packages
	Fitting in multi-dimensions - along energy axis, 3D etc. 
	XANES mapping -  alignment, tracking, PCA, cluster analysis
	Juptyper Notebook tools and examples
  
There are other python packages out there that perform these type of tasks -  
namely pymca, scikit-beam, mantis/spectromicroscopy and I would encourage you to look at them.

The focus of this package is on dealing with new experiments and problems and so it's based on using jupyter
and providing a toolkit to string together analysis rather than a smooth gui tool (similar to hyperspy)
	
flupy is a pure Python package, and so easy to install from source or with
``pip install flupy``.

I've tried to write a large number of examples to cover various use cases so that's probably the best place to start.

I've followed the excellent lmfit packages for this documentation and hyperspy/pymca to define a clear metadata structure.
Using the bug tracking software in GitHub Issues is encouraged for known problems and bug reports.
Please read `Contributing.md <.github/CONTRIBUTING.md>`_ before creating an Issue.


