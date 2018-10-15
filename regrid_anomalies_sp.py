import numpy
import hdf5storage

all_dat2 = numpy.fromfile('spah_daily365_1979-2017.dat', dtype='f8').reshape([14235, 88838])

matfiledata = {}
matfiledata['sp'] = all_dat2
hdf5storage.write(matfiledata, '.', 'data/sp01_365.mat', matlab_compatible=True)

all_dat5 = all_dat2.reshape([14235//5, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['sp'] = all_dat5
hdf5storage.write(matfiledata, '.', 'data/sp05_365.mat', matlab_compatible=True)

all_dat15 = all_dat5.reshape([2847//3, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['sp'] = all_dat15
hdf5storage.write(matfiledata, '.', 'data/sp15_365.mat', matlab_compatible=True)
