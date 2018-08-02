import numpy
import scipy.io
import xarray

variables = ['F', 'mspe', 'sig']
all_vars = []

for var in variables:
    for direction in ['sst_to_vort', 'vort_to_sst']:
        all_vars.append('{var}_{dir}'.format(var=var, dir=direction))


vars1 = dict(zip(all_vars, [(('lat', 'lon', 'delay'), numpy.full([241, 480, 31], numpy.NaN)) for _ in range(len(all_vars))]))
data1 = xarray.Dataset(vars1, coords={'lat': range(1, 242),
                                      'lon': range(1, 481),
                                      'delay': range(0, 31)})

vars5 = dict(zip(all_vars, [(('lat', 'lon', 'delay'), numpy.full((241, 480, 16), numpy.NaN)) for _ in range(len(all_vars))]))
data5 = xarray.Dataset(vars5, coords={'lat': range(1, 242), 'lon': range(1, 481), 'delay': range(0, 16)})

vars15 = dict(zip(all_vars, [(('lat', 'lon', 'delay'), numpy.full((241, 480, 11), numpy.NaN)) for _ in range(len(all_vars))]))
data15 = xarray.Dataset(vars15, coords={'lat': range(1, 242), 'lon': range(1, 481), 'delay': range(0, 11)})

for lon in range(1, 481):
    for var in variables:
        for direction in ['sst_to_vort', 'vort_to_sst']:
            name = '{var}_{dir}'.format(var=var, dir=direction)
            mat = scipy.io.loadmat('data/{name}_{lon}.mat'.format(name=name,
                                                                  lon=lon))[name][0]
            data1[name][:, lon - 1, :] = mat[0][0]
            data5[name][:, lon - 1, :] = mat[1][0]
            data15[name][:, lon - 1, :] = mat[2][0]

data1.to_netcdf('data1.nc')
data5.to_netcdf('data5.nc')
data15.to_netcdf('data15.nc')
