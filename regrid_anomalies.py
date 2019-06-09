import numpy
import hdf5storage

all_dat2 = numpy.fromfile('vortah850_daily365_1979-2017.dat',
                          dtype='f8').reshape([14235, 88838])

matfiledata = {}
matfiledata['vort'] = all_dat2
hdf5storage.write(matfiledata, '.', 'data/vort01_365.mat',
                  matlab_compatible=True)

all_dat5 = all_dat2.reshape([14235//5, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['vort'] = all_dat5
hdf5storage.write(matfiledata, '.', 'data/vort05_365.mat',
                  matlab_compatible=True)

all_dat2 = numpy.fromfile('sstah_daily365_1979-2017.dat',
                          dtype='f8').reshape([14235, 88838])

matfiledata = {}
matfiledata['sst'] = all_dat2
hdf5storage.write(matfiledata, '.', 'data/sst01_365.mat',
                  matlab_compatible=True)

all_dat5 = all_dat2.reshape([14235//5, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['sst'] = all_dat5
hdf5storage.write(matfiledata, '.', 'data/sst05_365.mat',
                  matlab_compatible=True)

for var in ['temp', 'q', 'sp', 'div']:
    all_dat = numpy.fromfile('{var}ah_daily365_1979-2017.dat'.format(var=var),
                             dtype='f8').reshape([14235, 88838])

    matfiledata = {}
    matfiledata[var] = all_dat
    hdf5storage.write(matfiledata, '.', 'data/{var}01_365.mat'.format(var),
                      matlab_compatible=True)

    all_dat.close()
