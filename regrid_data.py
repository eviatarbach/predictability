import pygrib
import numpy
import pandas

def leap_day(s):
    return (s.year % 4 == 0) & ((s.year % 100 != 0) | (s.year % 400 == 0)) & (s.month == 2) & (s.day == 29)

for var in ['div', 'sp', 'q', 'sst', 'temp', 'vort']:
    dat = pygrib.open('{var}.grib'.format(var=var)).select()

    rng = pandas.date_range(start='1/1/1979', end='12/31/2017', freq='D')

    mask = leap_day(rng)

    all_dat = numpy.zeros([56980, 88838])
    for i, dat2 in enumerate(dat):
        dat2.expand_grid(False)
        all_dat[i, :] = dat2.data()[0].data

    all_dat2 = all_dat.reshape([56980//4, -1, 88838]).mean(axis=1)
    all_dat2 = all_dat2[~mask, :]  # remove leap days
    all_dat2.tofile('{var}_daily365_1979-2017.dat'.format(var=var))

    dat.close()
