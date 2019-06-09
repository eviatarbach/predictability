import numpy

vort = numpy.fromfile('vortah850_daily365_1979-2017.dat',
                      dtype='f8').reshape([14235, 88838])
temp = numpy.fromfile('tempah_daily365_1979-2017.dat',
                      dtype='f8').reshape([14235, 88838])
q = numpy.fromfile('qah_daily365_1979-2017.dat',
                   dtype='f8').reshape([14235, 88838])
div = numpy.fromfile('divah_daily365_1979-2017.dat',
                     dtype='f8').reshape([14235, 88838])
sp = numpy.fromfile('spah_daily365_1979-2017.dat',
                    dtype='f8').reshape([14235, 88838])
sst = numpy.fromfile('sstah_daily365_1979-2017.dat',
                     dtype='f8').reshape([14235, 88838])

atmos_var = numpy.zeros(88838)

for cell in range(88838):
    X = numpy.zeros([5, 14235])
    X[0, :] = vort[:, cell]
    X[1, :] = temp[:, cell]
    X[2, :] = q[:, cell]
    X[3, :] = div[:, cell]
    X[4, :] = sp[:, cell]
    atmos_var[cell] = numpy.linalg.det(numpy.cov(X))

numpy.save('var_atmos.npy', atmos_var)

sst_var = numpy.var(sst, axis=0)
numpy.save('var_sst.npy', sst_var)
