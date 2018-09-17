import numpy
import hdf5storage

all_dat2 = numpy.fromfile('data/vort01_365.dat', dtype='f8').reshape([14235, 88838])
matfiledata = {}
matfiledata['vort'] = all_dat2
hdf5storage.write(matfiledata, '.', 'data/vort01_365.mat', matlab_compatible=True)

all_dat2 = numpy.fromfile('data/sst01_365.dat', dtype='f8').reshape([14235, 88838])
matfiledata = {}
matfiledata['sst'] = all_dat2
hdf5storage.write(matfiledata, '.', 'data/sst01_365.mat', matlab_compatible=True)
