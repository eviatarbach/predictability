import numpy
import scipy.io
import xarray

variables = ['F', 'mspe_full', 'mspe_reduced', 'sig']
all_vars = ['times']

for var in variables:
    for direction in ['sst_to_atmos', 'atmos_to_sst']:
        all_vars.append('{var}_{dir}'.format(var=var, dir=direction))


vars1 = dict(zip(all_vars, [(('cell', 'delay'), numpy.full((88838, 361), numpy.NaN)) for _ in range(len(all_vars))]))
data1 = xarray.Dataset(vars1, coords={'cell': range(1, 88839),
                                      'delay': range(0, 361)})

for offset in range(1, 88838, 60):
    end = min(88838 - offset + 1, 60)
    for var in variables:
        for direction in ['sst_to_atmos', 'atmos_to_sst']:
            name = '{var}_{dir}'.format(var=var, dir=direction)
            try:
                mat = scipy.io.loadmat('data_atmos/{name}_{cell}.mat'.format(name=name,
                                                                       cell=offset))[name][0]
                data1[name][offset - 1:offset + 60 - 1, :] = mat[0][:end, :]
            except Exception as a:
                pass
                #print(offset)
                #print(a)
    try:
        mat = scipy.io.loadmat('data_atmos/times_atmos_{cell}.mat'.format(cell=offset))['times'][0]
        data1['times'][offset - 1:offset + 60 - 1, 0] = mat[0][:end, 0]
    except Exception as a:
        print(offset)

data1.to_netcdf('data_atmos/data01.nc')
