import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import cmocean
import scipy.io
import numpy

F_sst_to_vort = scipy.io.loadmat('data/F_sst_to_vort.mat')['F_sst_to_vort'][0][1]
F_vort_to_sst = scipy.io.loadmat('data/F_vort_to_sst.mat')['F_vort_to_sst'][0][1]

for var in ['ratio', 'atmos', 'ocean', 'both']:
    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, color='#D9D9D9')
    ax.add_feature(cfeature.COASTLINE, linewidth=0.3)

    lat = numpy.arange(-90, -90 + 241*0.75, 0.75)
    lon = numpy.arange(0, 480*0.75, 0.75)

    if var == 'ratio':
        ratios = numpy.log(F_vort_to_sst/F_sst_to_vort)
        ratios[ratios > 5] = 5
        ratios[ratios < -5] = -5
        var_to_plot = ratios
        levels = numpy.linspace(-5, 5, 40)
        cmap = cmocean.cm.balance
    elif var == 'atmos':
        var_to_plot = F_vort_to_sst
        levels = numpy.linspace(0, numpy.nanmax(F_vort_to_sst), 40)
        cmap = cmocean.cm.amp
    elif var == 'ocean':
        var_to_plot = F_sst_to_vort
        levels = numpy.linspace(0, numpy.nanmax(F_sst_to_vort), 40)
        upper_half = numpy.linspace(0.5, 1, 40)
        cmap = cmocean.cm.balance_r
        colors = cmap(upper_half)
        cmap = LinearSegmentedColormap.from_list('Upper Half', colors)
    elif var == 'both':
        both = F_vort_to_sst.copy()
        both[F_sst_to_vort > F_vort_to_sst] = -F_sst_to_vort[F_sst_to_vort > F_vort_to_sst] 
        var_to_plot = both
        levels = numpy.linspace(-0.15, 0.15, 40)
        cmap = cmocean.cm.balance

    var_cyc, lon_cyc = add_cyclic_point(var_to_plot.T, coord=lon)
    
    latt, lonn = numpy.meshgrid(lat, lon_cyc)

    plt.contourf(lonn, latt, var_cyc.T, vmin=min(levels), vmax=max(levels),
                 cmap=cmap,
                 levels=levels,
                 transform=ccrs.PlateCarree())
    plt.colorbar()

    plt.savefig('map_comparison_{var}.eps'.format(var=var))

    plt.gcf().clear()
