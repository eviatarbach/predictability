import pygrib
import numpy
import hdf5storage

dat = pygrib.open('/lustre/ebach/sfc.grib').select()

all_dat = numpy.zeros([56980, 88838])
for i, dat2 in enumerate(dat):
    dat2.expand_grid(False)
    all_dat[i, :] = dat2.data()[0].data

all_dat2 = all_dat.reshape([56980//4, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['sst'] = all_dat2
hdf5storage.write(matfiledata, '.', 'data/sst01_g.mat', matlab_compatible=True)

all_dat5 = all_dat2.reshape([14245//5, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['sst'] = all_dat5
hdf5storage.write(matfiledata, '.', 'data/sst05_g.mat', matlab_compatible=True)

all_dat15 = all_dat5[:2847, :].reshape([2847//3, -1, 88838]).mean(axis=1)
matfiledata = {}
matfiledata['sst'] = all_dat15
hdf5storage.write(matfiledata, '.', 'data/sst15_g.mat', matlab_compatible=True)
