import numpy
import scipy.io
import xarray

variables = ['F', 'sig']
all_vars = []

for var in variables:
    for direction in ['sst_to_vort', 'vort_to_sst']:
        all_vars.append('{var}_{dir}'.format(var=var, dir=direction))

vars_winter = dict(zip(all_vars, [(('cell'), numpy.full(88838, numpy.NaN)) for _ in range(len(all_vars))]))
data_winter = xarray.Dataset(vars_winter, coords={'cell': range(1, 88839)})

vars_spring = dict(zip(all_vars, [(('cell'), numpy.full(88838, numpy.NaN)) for _ in range(len(all_vars))]))
data_spring = xarray.Dataset(vars_spring, coords={'cell': range(1, 88839)})

vars_summer = dict(zip(all_vars, [(('cell'), numpy.full(88838, numpy.NaN)) for _ in range(len(all_vars))]))
data_summer = xarray.Dataset(vars_summer, coords={'cell': range(1, 88839)})

vars_fall = dict(zip(all_vars, [(('cell'), numpy.full(88838, numpy.NaN)) for _ in range(len(all_vars))]))
data_fall = xarray.Dataset(vars_fall, coords={'cell': range(1, 88839)})

for offset in range(1, 88838, 185):
    end = min(88838 - offset + 1, 185)
    for var in variables:
        for direction in ['sst_to_vort', 'vort_to_sst']:
            name = '{var}_seasonal_{dir}'.format(var=var, dir=direction)
            name2 = '{var}_{dir}'.format(var=var, dir=direction)
            try:
                mat = scipy.io.loadmat('data/{name}_{cell}.mat'.format(name=name,
                                                                       cell=offset))[name2][0]
                data_winter[name2][offset - 1:offset + 185 - 1] = mat[0][:, 0][:end]
                data_spring[name2][offset - 1:offset + 185 - 1] = mat[1][:, 0][:end]
                data_summer[name2][offset - 1:offset + 185 - 1] = mat[2][:, 0][:end]
                data_fall[name2][offset - 1:offset + 185 - 1] = mat[3][:, 0][:end]
            except Exception as a:
                print(offset)
                print(a)

data_winter.to_netcdf('data/data_winter.nc')
data_spring.to_netcdf('data/data_spring.nc')
data_summer.to_netcdf('data/data_summer.nc')
data_fall.to_netcdf('data/data_fall.nc')
