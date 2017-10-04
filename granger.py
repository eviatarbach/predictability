import multiprocessing

import numpy
from statsmodels.tsa.stattools import grangercausalitytests

sst_daily = numpy.fromfile('data/sstah_eraint_0.75x0.75_daily365_1979-2016.dat', numpy.float32)
sst_daily = sst_daily.reshape([13870, 241, 480])

vort_daily = numpy.fromfile('data/vort850ahocn_eraint_0.75x0.75_daily365_1979-2016.dat',
                            numpy.float32)
vort_daily = vort_daily.reshape([13870, 241, 480])


def calc_ratio(ts):
    if ts[0, 0] == -9999:
        ratio = 0
    else:
        F_1, p_1, _ = grangercausalitytests(ts, maxlag=10, verbose=False)[10][0]['lrtest']
        F_2, p_2, _ = grangercausalitytests(numpy.flipud(ts), maxlag=10, verbose=False)[10][0]['lrtest']
        if (p_1 < 0.1) or (p_2 < 0.1):
            if F_1 > F_2:
                ratio = min(F_1/F_2, 15)
            else:
                ratio = max(-F_2/F_1, -15)
        else:
            ratio = 0
    print(ratio)
    return ratio


# Lazy way to get data
idx = numpy.ndindex(241, 480)
data_getter = map(lambda i: numpy.vstack([sst_daily[:, i[0], i[1]], vort_daily[:, i[0], i[1]]]).T, idx)

pool = multiprocessing.Pool(processes=20)
data = pool.map(calc_ratio, data_getter, chunksize=36)
