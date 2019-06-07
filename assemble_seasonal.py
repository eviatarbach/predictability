import numpy
import scipy.io
import xarray

variables = ['F', 'mspe_full', 'mspe_reduced', 'sig']
all_vars = []

for var in variables:
    for direction in ['sst_to_atmos', 'atmos_to_sst']:
        all_vars.append('{var}_{dir}'.format(var=var, dir=direction))

vars_winter = dict(zip(all_vars, [(('cell'), numpy.full(88838, numpy.NaN)) for _ in range(len(all_vars))]))
data_winter = xarray.Dataset(vars_winter, coords={'cell': range(1, 88839)})

vars_spring = dict(zip(all_vars, [(('cell'), numpy.full(88838, numpy.NaN)) for _ in range(len(all_vars))]))
data_spring = xarray.Dataset(vars_spring, coords={'cell': range(1, 88839)})

vars_summer = dict(zip(all_vars, [(('cell'), numpy.full(88838, numpy.NaN)) for _ in range(len(all_vars))]))
data_summer = xarray.Dataset(vars_summer, coords={'cell': range(1, 88839)})

vars_fall = dict(zip(all_vars, [(('cell'), numpy.full(88838, numpy.NaN)) for _ in range(len(all_vars))]))
data_fall = xarray.Dataset(vars_fall, coords={'cell': range(1, 88839)})

missing = []
for offset in range(1, 88838, 60):
    end = min(88838 - offset + 1, 60)
    for var in variables:
        for direction in ['sst_to_atmos', 'atmos_to_sst']:
            name = '{var}_seasonal_{dir}'.format(var=var, dir=direction)
            name2 = '{var}_{dir}'.format(var=var, dir=direction)
            try:
                mat = scipy.io.loadmat('data_atmos/{name}_{cell}.mat'.format(name=name,
                                                                       cell=offset))[name2][0]
                data_winter[name2][offset - 1:offset + 60 - 1] = mat[0][:, 0][:end]
                data_spring[name2][offset - 1:offset + 60 - 1] = mat[1][:, 0][:end]
                data_summer[name2][offset - 1:offset + 60 - 1] = mat[2][:, 0][:end]
                data_fall[name2][offset - 1:offset + 60 - 1] = mat[3][:, 0][:end]
            except Exception as a:
                if name2 == 'F_sst_to_atmos':
                    missing.append(offset)
                    print(offset)

data_winter.to_netcdf('data_atmos/data_winter.nc')
data_spring.to_netcdf('data_atmos/data_spring.nc')
data_summer.to_netcdf('data_atmos/data_summer.nc')
data_fall.to_netcdf('data_atmos/data_fall.nc')
