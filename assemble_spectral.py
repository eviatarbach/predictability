import numpy
import scipy.io
import xarray

variables = ['f']
all_vars = []

for var in variables:
    for direction in ['sst_to_atmos', 'atmos_to_sst']:
        all_vars.append('{var}_{dir}'.format(var=var, dir=direction))


vars = dict(zip(all_vars, [(('cell', 'freq'), numpy.full((88838, 301), numpy.NaN)) for _ in range(len(all_vars))]))
data = xarray.Dataset(vars, coords={'cell': range(1, 88839),
                                    'freq': range(301)})

for offset in range(1, 88838, 60):
    end = min(88838 - offset + 1, 60)
    for var in variables:
        for direction in ['sst_to_atmos', 'atmos_to_sst']:
            name = '{var}_{dir}'.format(var=var, dir=direction)
            try:
                mat = scipy.io.loadmat('data_atmos/{name}_{cell}.mat'.format(name=name,
                                                                       cell=offset))[name]
                data[name][offset - 1:offset + 60 - 1, :] = mat[:end, :]
            except Exception as a:
                pass
                print(offset)

data.to_netcdf('data_atmos/spectral.nc')
