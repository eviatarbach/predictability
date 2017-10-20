import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import scipy.io
import numpy

ax = plt.axes(projection=ccrs.Robinson())
ax.coastlines()

lat = numpy.arange(-90, -90 + 241*0.75, 0.75)
lon = numpy.arange(0, 480*0.75, 0.75)

latt, lonn = numpy.meshgrid(lat, lon)

#ratios = numpy.load('data/granger.npy').reshape([241, 480]).T
ratios = scipy.io.loadmat('granger.mat')['ratios']

ratios[ratios == 0] = numpy.nan

plt.contourf(lonn, latt, ratios, vmin=-15, vmax=15, cmap='seismic', levels=numpy.linspace(-15, 15, 10)[1:-1])
