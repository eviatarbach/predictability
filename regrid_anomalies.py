import numpy
import hdf5storage

or var in ['vort', 'sst', 'temp', 'q', 'sp', 'div']:
    all_dat = numpy.fromfile('{var}ah_daily365_1979-2017.dat'.format(var=var),
                             dtype='f8').reshape([14235, 88838])

    matfiledata = {}
    matfiledata[var] = all_dat
    hdf5storage.write(matfiledata, '.', 'data/{var}01_365.mat'.format(var),
                      matlab_compatible=True)

    all_dat.close()
