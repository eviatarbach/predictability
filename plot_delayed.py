import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import scipy.io
import numpy

sig_90_sst_to_vort = scipy.io.loadmat('data/sig_90_sst_to_vort_delay.mat')['sig_90_sst_to_vort'][0]
sig_90_vort_to_sst = scipy.io.loadmat('data/sig_90_vort_to_sst_delay.mat')['sig_90_vort_to_sst'][0]
sig_95_sst_to_vort = scipy.io.loadmat('data/sig_95_sst_to_vort_delay.mat')['sig_95_sst_to_vort'][0]
sig_95_vort_to_sst = scipy.io.loadmat('data/sig_95_vort_to_sst_delay.mat')['sig_95_vort_to_sst'][0]
F_sst_to_vort = scipy.io.loadmat('data/F_sst_to_vort_delay.mat')['F_sst_to_vort'][0]
F_vort_to_sst = scipy.io.loadmat('data/F_vort_to_sst_delay.mat')['F_vort_to_sst'][0]

for i in range(3):
    for j in range(len(F_vort_to_sst[i][0])):
        plt.figure(figsize=(6, 3))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()

        lat = numpy.arange(-90, -90 + 241*0.75, 0.75)
        lon = numpy.arange(0, 480*0.75, 0.75)

        latt, lonn = numpy.meshgrid(lat, lon)

        ratios = F_vort_to_sst[i][0][j]/F_sst_to_vort[i][0][j]
        ratios[ratios < 1] = -1/ratios[ratios < 1]

        sig_90_vort_to_sst[i][0][j][numpy.isnan(sig_90_vort_to_sst[i][0][j])] = 0
        sig_90_sst_to_vort[i][0][j][numpy.isnan(sig_90_sst_to_vort[i][0][j])] = 0

        sig_95_vort_to_sst[i][0][j][numpy.isnan(sig_95_vort_to_sst[i][0][j])] = 0
        sig_95_sst_to_vort[i][0][j][numpy.isnan(sig_95_sst_to_vort[i][0][j])] = 0

        sig = sig_95_sst_to_vort[i][0][j].astype(bool) | sig_95_vort_to_sst[i][0][j].astype(bool)
        ratios[~sig] = numpy.nan

        plt.contourf(lonn, latt, ratios, vmin=-15, vmax=15, cmap='seismic', levels=numpy.hstack([numpy.linspace(-15, -1, 4), numpy.linspace(1, 15, 4)]), transform=ccrs.PlateCarree())
        plt.colorbar()

        plt.savefig('map_delay_{i}_{j}.eps'.format(i=i, j=j))

        plt.gcf().clear()
