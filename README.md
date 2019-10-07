# Code for "Local atmosphere–ocean predictability: dynamical origins, lead times, and seasonality"

This repository contains the code for the paper "Local atmosphere–ocean predictability: dynamical origins, lead times, and seasonality" by Eviatar Bach, Safa Motesharrei, Eugenia Kalnay, and Alfredo Ruiz-Barradas.

All the code was written by Eviatar Bach, except for the anomaly time-series code which was written by Alfredo Ruiz-Barradas. All code by Eviatar Bach is licensed under the GNU Public License v3.0.

The Granger causality analysis is done using the MVGC Matlab library, and plotting is done in Python.

You can contact me with any questions at eviatarbach@protonmail.com.

## Output of analysis
The model output for the main analysis can be downloaded [here](http://www.terpconnect.umd.edu/~ebach/data01.nc), and for the spectral [here](http://www.terpconnect.umd.edu/~ebach/spectral.nc).

## Libraries used

Matlab:
- [MVGC](http://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html)

Python:
- [cartopy](http://scitools.org.uk/cartopy/)
- [cmocean](https://matplotlib.org/cmocean/)
- [ecmwfapi](https://pypi.org/project/ecmwf-api-client/)
- [hdf5storage](https://pythonhosted.org/hdf5storage/)
- [Matplotlib](https://matplotlib.org/)
- [NumPy](http://www.numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [pygrib](https://github.com/jswhit/pygrib)
- [pyresample](https://pyresample.readthedocs.io/en/latest/)
- [SciPy](https://scipy.org/scipylib/index.html)
- [statsmodels](https://www.statsmodels.org/stable/index.html)
- [xarray](http://xarray.pydata.org/en/stable/)

You also need [Jupyter](https://jupyter.org/) for the notebook that generates the plots.

## Description of files

The files are listed in the order that they should be run:

- **daily_data.py**: Takes daily means of the reanalysis data and removes leap days.
- **retrieve_data.py**: Retrieves the reanalysis data (ERA-Interim) from ECMWF.
- **anomalies_vort.f**: Fortran program to remove the seasonal and subseasonal components from the vorticity data, leaving the anomalies. Replace vorticity for the other variables.
- **convert_anomalies.py**: Convert the anomalies to HDF5 for Matlab.
- **granger.m**: Performs the main Granger causality analysis, as well as that for lead time (generates output used in Figs. 2–5 and 8a).
- **granger_seasonal.m**: Performs the Granger causality analysis by season (generates output used in Fig. 7).
- **granger_spectral.m**: Performs the spectral Granger causality analysis (generates output used in Figs. 6, 8b, and S3 in the Supplementary Information)
- **assemble.py**: Assembles the pieces after ``granger.m`` into a NetCDF file.
- **assemble_spectral.py**: Assembles the pieces after ``granger_spectral.m`` into a NetCDF file.
- **assemble_seasonal.py**: Assembles the pieces after ``granger_seasonal.m`` into NetCDF files.
- **variance.py**: Computes the variance of the SST data as well as the generalized variance of the atmospheric variables.
- **Plots.ipynb**: Jupyter notebook that generates all the figures in the paper (except Figs. 1, 9, and S1 and S2 in the Supplementary Information).
