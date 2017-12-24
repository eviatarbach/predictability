import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
from matplotlib import ticker
import cmocean
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
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, color='#D9D9D9')
    ax.add_feature(cfeature.COASTLINE, linewidth=0.3)

    lat = numpy.arange(-90, -90 + 241*0.75, 0.75)
    lon = numpy.arange(0, 480*0.75, 0.75)

    ratios = numpy.log(F_vort_to_sst[i]/F_sst_to_vort[i])

    sig_90_vort_to_sst[i][numpy.isnan(sig_90_vort_to_sst[i])] = 0
    sig_90_sst_to_vort[i][numpy.isnan(sig_90_sst_to_vort[i])] = 0

    sig_95_vort_to_sst[i][numpy.isnan(sig_95_vort_to_sst[i])] = 0
    sig_95_sst_to_vort[i][numpy.isnan(sig_95_sst_to_vort[i])] = 0

    sig = sig_95_sst_to_vort[i].astype(bool) | sig_95_vort_to_sst[i].astype(bool)
    ratios[~sig] = numpy.nan

    ratios_cyc, lon_cyc = add_cyclic_point(ratios.T, coord=lon)
    
    latt, lonn = numpy.meshgrid(lat, lon_cyc)

    plt.contourf(lonn, latt, ratios_cyc.T, vmin=-5, vmax=5,
                 cmap=cmocean.cm.balance, levels=numpy.linspace(-5, 5, 40),
                 transform=ccrs.PlateCarree(), extend='both')
    cb = plt.colorbar(orientation='horizontal', fraction=0.05, pad=0.04)
    cb.set_label(r'$\alpha$')
    tick_locator = ticker.MaxNLocator(nbins=9)
    cb.locator = tick_locator
    cb.update_ticks()

    plt.savefig('map_{i}.eps'.format(i=i))

    plt.gcf().clear()
