import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import ticker
import cmocean
import scipy.io
import numpy

xlabels = [r'$\alpha$', r'$\mathrm{GC}_{\mathrm{Vort}\rightarrow\mathrm{SST}}$',
           r'$\mathrm{GC}_{\mathrm{SST}\rightarrow\mathrm{Vort}}$',
           r'$\operatorname{sgn}(\mathrm{GC}_{\mathrm{Vort}\rightarrow\mathrm{SST}}-\mathrm{GC}_{\mathrm{SST}\rightarrow\mathrm{Vort}})\,\max(\mathrm{GC}_{\mathrm{SST}\rightarrow\mathrm{Vort}}, \mathrm{GC}_{\mathrm{Vort}\rightarrow\mathrm{SST}})$']

F_sst_to_vort = scipy.io.loadmat('data/F_sst_to_vort.mat')['F_sst_to_vort'][0]
F_vort_to_sst = scipy.io.loadmat('data/F_vort_to_sst.mat')['F_vort_to_sst'][0]

for j, interval in enumerate(['daily', 'pentad', '15_day']):
    F_sst_to_vort_j = F_sst_to_vort[j]
    F_vort_to_sst_j = F_vort_to_sst[j]
    for i, var in enumerate(['ratio', 'atmos', 'ocean', 'both']):
        plt.figure()
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, color='#D9D9D9')
        ax.add_feature(cfeature.COASTLINE, linewidth=0.3)

        lat = numpy.arange(-90, -90 + 241*0.75, 0.75)
        lon = numpy.arange(0, 480*0.75, 0.75)

        if var == 'ratio':
            ratios = numpy.log(F_vort_to_sst_j/F_sst_to_vort_j)
            var_to_plot = ratios
            levels = numpy.linspace(-5, 5, 40)
            cmap = cmocean.cm.balance
        elif var == 'atmos':
            var_to_plot = F_vort_to_sst_j
            if interval == 'pentad':
                levels = numpy.linspace(0, 0.15, 40)
            else:
                levels = numpy.linspace(0, 0.05, 40)
            cmap = cmocean.cm.amp
        elif var == 'ocean':
            var_to_plot = F_sst_to_vort_j
            if interval == 'pentad':
                levels = numpy.linspace(0, 0.04, 40)
            else:
                levels = numpy.linspace(0, 0.01, 40)
            upper_half = numpy.linspace(0.5, 1, 40)
            cmap = cmocean.cm.balance_r
            colors = cmap(upper_half)
            cmap = LinearSegmentedColormap.from_list('Upper Half', colors)
        elif var == 'both':
            both = F_vort_to_sst_j.copy()
            both[F_sst_to_vort_j > F_vort_to_sst_j] = -F_sst_to_vort_j[F_sst_to_vort_j > F_vort_to_sst_j] 
            var_to_plot = both
            levels = numpy.linspace(-0.15, 0.15, 40)
            cmap = cmocean.cm.balance

        var_cyc, lon_cyc = add_cyclic_point(var_to_plot.T, coord=lon)
        
        latt, lonn = numpy.meshgrid(lat, lon_cyc)

        plt.contourf(lonn, latt, var_cyc.T, vmin=min(levels), vmax=max(levels),
                     cmap=cmap, levels=levels, transform=ccrs.PlateCarree(),
                     extend='both')
        cb = plt.colorbar(orientation='horizontal', fraction=0.05, pad=0.04)
        cb.set_label(xlabels[i])
        tick_locator = ticker.MaxNLocator(nbins=9)
        cb.locator = tick_locator
        cb.update_ticks()

        #plt.title('{var} forcing using Granger causality ({resolution} resolution)'.format(var=['Oceanic', 'Atmospheric'][j], resolution=['daily', 'pentad'][i]))
        plt.savefig('map_comparison_{var}_{interval}.eps'.format(var=var, interval=interval))

        plt.gcf().clear()
