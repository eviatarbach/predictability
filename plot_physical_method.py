import numpy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import ticker
import cmocean

daily = numpy.fromfile('data/freq_sstvort850ah_eraint_0.75x0.75_daily365_1979-2016.dat', numpy.float32)
daily[daily == -9999] = numpy.nan
daily = daily.reshape([3, 4, 241, 480])

pentad = numpy.fromfile('data/freq_sstvort850ah_eraint_0.75x0.75_pentad_1979-2016.dat', numpy.float32)
pentad[pentad == -9999] = numpy.nan
pentad = pentad.reshape([3, 4, 241, 480])

for i, data in enumerate([daily, pentad]):
    for j, var in enumerate(['ocean', 'atmos', 'ratio']):
        ocean = data[0, 0, :, :] + data[0, 1, :, :]
        atmos = data[0, 2, :, :] + data[0, 3, :, :]

        plt.figure()
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, color='#D9D9D9')
        ax.add_feature(cfeature.COASTLINE, linewidth=0.3)

        lat = numpy.arange(-90, -90 + 241*0.75, 0.75)
        lon = numpy.arange(0, 480*0.75, 0.75)

        if var == 'atmos':
            var_to_plot = atmos
            levels = numpy.linspace(numpy.nanmedian(var_to_plot),
                                    numpy.nanmax(var_to_plot), 40)
            cmap = cmocean.cm.amp
        elif var == 'ocean':
            var_to_plot = ocean
            levels = numpy.linspace(numpy.nanmedian(var_to_plot),
                                    numpy.nanmax(var_to_plot), 40)
            upper_half = numpy.linspace(0.5, 1, 40)
            cmap = cmocean.cm.balance_r
            colors = cmap(upper_half)
            cmap = LinearSegmentedColormap.from_list('Upper Half', colors)
        elif var == 'ratio':
            var_to_plot = numpy.log(atmos/ocean)
            levels = numpy.linspace(-1, 1, 40)
            cmap = cmocean.cm.balance

        var_cyc, lon_cyc = add_cyclic_point(var_to_plot, coord=lon)

        latt, lonn = numpy.meshgrid(lat, lon_cyc)

        plt.contourf(lonn, latt, var_cyc.T, vmin=min(levels),
                     vmax=max(levels), cmap=cmap, levels=levels,
                     transform=ccrs.PlateCarree(), extend='both')
        cb = plt.colorbar(orientation='horizontal', fraction=0.05, pad=0.04)
        cb.set_label('Anomaly count')
        tick_locator = ticker.MaxNLocator(nbins=9)
        cb.locator = tick_locator
        cb.update_ticks()

        plt.title('{var} forcing using dynamical rule ({resolution} resolution)'.format(var=['Oceanic', 'Atmospheric', 'Ratio'][j], resolution=['daily', 'pentad'][i]))
        plt.savefig('map_physical_{resolution}_{var}.eps'.format(resolution=['daily', 'pentad'][i], var=var))

        plt.gcf().clear()
