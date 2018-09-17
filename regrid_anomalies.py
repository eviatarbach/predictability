import numpy
import hdf5storage

all_dat2 = numpy.fromfile('vortah850_daily365_1979-2017.dat', dtype='f8').reshape([14235, 88838])

matfiledata = {}
matfiledata['vort'] = all_dat2
hdf5storage.write(matfiledata, '.', 'data/vort01_365.mat', matlab_compatible=True)

all_dat5 = all_dat2.reshape([14235//5, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['vort'] = all_dat5
hdf5storage.write(matfiledata, '.', 'data/vort05_365.mat', matlab_compatible=True)

all_dat15 = all_dat5.reshape([2847//3, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['vort'] = all_dat15
hdf5storage.write(matfiledata, '.', 'data/vort15_365.mat', matlab_compatible=True)

all_dat2 = numpy.fromfile('sstah_daily365_1979-2017.dat', dtype='f8').reshape([14235, 88838])

matfiledata = {}
matfiledata['sst'] = all_dat2
hdf5storage.write(matfiledata, '.', 'data/sst01_365.mat', matlab_compatible=True)

all_dat5 = all_dat2.reshape([14235//5, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['sst'] = all_dat5
hdf5storage.write(matfiledata, '.', 'data/sst05_365.mat', matlab_compatible=True)

all_dat15 = all_dat5.reshape([2847//3, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['sst'] = all_dat15
hdf5storage.write(matfiledata, '.', 'data/sst15_365.mat', matlab_compatible=True)