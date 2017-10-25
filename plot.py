import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import scipy.io
import numpy

sig_90_sst_to_vort = scipy.io.loadmat('data/sig_90_sst_to_vort.mat')['sig_90_sst_to_vort'][0]
sig_90_vort_to_sst = scipy.io.loadmat('data/sig_90_vort_to_sst.mat')['sig_90_vort_to_sst'][0]
sig_95_sst_to_vort = scipy.io.loadmat('data/sig_95_sst_to_vort.mat')['sig_95_sst_to_vort'][0]
sig_95_vort_to_sst = scipy.io.loadmat('data/sig_95_vort_to_sst.mat')['sig_95_vort_to_sst'][0]
F_sst_to_vort = scipy.io.loadmat('data/F_sst_to_vort.mat')['F_sst_to_vort'][0]
F_vort_to_sst = scipy.io.loadmat('data/F_vort_to_sst.mat')['F_vort_to_sst'][0]
times = scipy.io.loadmat('data/times.mat')['times'][0]

for i in range(3):
	plt.figure(figsize=(6, 3))
	ax = plt.axes(projection=ccrs.PlateCarree())
	ax.coastlines()

	lat = numpy.arange(-90, -90 + 241*0.75, 0.75)
	lon = numpy.arange(0, 480*0.75, 0.75)

	latt, lonn = numpy.meshgrid(lat, lon)

	ratios = F_vort_to_sst[i]/F_sst_to_vort[i]
	ratios[ratios < 1] = -1/ratios[ratios < 1]
	
	sig_90_vort_to_sst[i][numpy.isnan(sig_90_vort_to_sst[i])] = 0
	sig_90_sst_to_vort[i][numpy.isnan(sig_90_sst_to_vort[i])] = 0

	sig_95_vort_to_sst[i][numpy.isnan(sig_95_vort_to_sst[i])] = 0
	sig_95_sst_to_vort[i][numpy.isnan(sig_95_sst_to_vort[i])] = 0

	sig = sig_95_sst_to_vort[i].astype(bool) | sig_95_vort_to_sst[i].astype(bool)
	ratios[~sig] = numpy.nan

	plt.contourf(lonn, latt, ratios, vmin=-15, vmax=15, cmap='seismic', levels=numpy.linspace(-15, 15, 10)[1:-1], transform=ccrs.PlateCarree())
	plt.colorbar()

	plt.savefig('map_{i}.eps'.format(i=i))

	plt.gcf().clear()
