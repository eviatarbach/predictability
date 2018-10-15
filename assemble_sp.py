import numpy
import scipy.io
import xarray

variables = ['F', 'mspe', 'sig']
all_vars = ['times']

for var in variables:
    for direction in ['sst_to_sp', 'sp_to_sst']:
        all_vars.append('{var}_{dir}'.format(var=var, dir=direction))


vars1 = dict(zip(all_vars, [(('cell', 'delay'), numpy.full((88838, 31), numpy.NaN)) for _ in range(len(all_vars))]))
data1 = xarray.Dataset(vars1, coords={'cell': range(1, 88839),
                                      'delay': range(0, 31)})

vars5 = dict(zip(all_vars, [(('cell', 'delay'), numpy.full((88838, 16), numpy.NaN)) for _ in range(len(all_vars))]))
data5 = xarray.Dataset(vars5, coords={'cell': range(1, 88839), 'delay': range(0, 16)})

vars15 = dict(zip(all_vars, [(('cell', 'delay'), numpy.full((88838, 11), numpy.NaN)) for _ in range(len(all_vars))]))
data15 = xarray.Dataset(vars15, coords={'cell': range(1, 88839), 'delay': range(0, 11)})

for offset in range(1, 88838, 185):
    end = min(88838 - offset + 1, 185)
    for var in variables:
        for direction in ['sst_to_sp', 'sp_to_sst']:
            name = '{var}_{dir}'.format(var=var, dir=direction)
            try:
                mat = scipy.io.loadmat('data_sp/{name}_{cell}.mat'.format(name=name,
                                                                       cell=offset))[name][0]
                data1[name][offset - 1:offset + 185 - 1, :] = mat[0][:end, :]
                data5[name][offset - 1:offset + 185 - 1, :] = mat[1][:end, :]
                data15[name][offset - 1:offset + 185 - 1, :] = mat[2][:end, :]
            except Exception as a:
                pass
                #print(offset)
                #print(a)
    try:
        mat = scipy.io.loadmat('data_sp/times_sp_{cell}.mat'.format(cell=offset))['times'][0]
        data1['times'][offset - 1:offset + 185 - 1, 0] = mat[0][:end, 0]
        data5['times'][offset - 1:offset + 185 - 1, 0] = mat[1][:end, 0]
        data15['times'][offset - 1:offset + 185 - 1, 0] = mat[2][:end, 0]
    except Exception as a:
        print(offset)

data1.to_netcdf('data_sp/data01.nc')
data5.to_netcdf('data_sp/data05.nc')
data15.to_netcdf('data_sp/data15.nc')
