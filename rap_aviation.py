#!/usr/bin/python3

import pygrib
import metpy.calc as mpcalc
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from scipy.ndimage.filters import gaussian_filter as gfilter

# function for computing the magnitude of a vector field
def windSpeed(u,v):
    return numpy.sqrt(u**2 + v**2)

# function for performing contour plots
def contourPlot(var1,var2,lev):
    assert (var1.shape == var2.shape),"Variables must be the same shape!"
    # create a new figure
    plt.clf()
    fig = plt.gcf()
    fig.set_size_inches(12.5,10.0)
    # set up the basemap and grid
    m = plotBasemap()
    ny,nx = var1.shape[0],var1.shape[1]
    longitudes,latitudes = m.makegrid(nx,ny)
    x,y = m(longitudes,latitudes)
    # create the filled contour plot
    cs_fill = m.contourf(x,y,var1,cmap='gist_ncar',levels=numpy.arange(30,150,5))
    cs_line = m.contour(x,y,var2,levels=numpy.arange(4,40,2),colors='k',linewidths=0.5)
    cbar = m.colorbar(cs_fill,location='bottom',pad='5%')
    cbar.set_label('Wind Speed (kt)')
    plt.clabel(cs_line,cs_line.levels,inline=True,fmt='%.0f',fontsize=8)
    plt.title('RAP Analyzed %d mb Wind Speed and Ellrod Index' % lev)
    plt.savefig('plots/rap_wspd_%.0f.png' % lev,bbox_inches='tight')

# function for creating a basemap
def plotBasemap():
    m = Basemap(llcrnrlon=-126.138,llcrnrlat=16.281,urcrnrlon=-57.383,urcrnrlat=55.481,projection='lcc',lat_1=25.0,lat_2=25.0,lon_0=-95.0,resolution ='i')
    m.drawcoastlines(color='grey')
    m.drawstates(color='grey')
    m.drawcountries(color='grey')
    return m

# function for computing deformation and convergence
def windProducts(x,y,u,v):
    dx,dy = mpcalc.lat_lon_grid_spacing(x,y)
    deformation = numpy.zeros(u.shape)
    convergence = numpy.zeros(u.shape)
    for i in range(u.shape[0]):
        shear,stretch = mpcalc.shearing_stretching_deformation(u[i,:,:],v[i,:,:],dx,dy)
        deformation[i] = numpy.sqrt(shear**2 + stretch**2)
        convergence[i] = -mpcalc.divergence(u[i,:,:],v[i,:,:],dx,dy)
    return deformation,convergence

# function for approximating vertical wind shear
def windShear(u,v,z):
    vws = numpy.zeros(u.shape)
    for i in range(u.shape[0]):
        if i > 0 and i < u.shape[0]-1:
            # first we compute the wind shear with the higher level
            du21 = u[i-1,:,:] - u[i,:,:]
            dv21 = v[i-1,:,:] - v[i,:,:]
            dz21 = z[i-1,:,:] - z[i,:,:]
            vws21 = numpy.sqrt(du21**2 + dv21**2) / dz21
            # then we compute the wind shear with the lower level
            du32 = u[i,:,:] - u[i+1,:,:]
            dv32 = v[i,:,:] - v[i+1,:,:]
            dz32 = z[i,:,:] - z[i+1,:,:]
            vws32 = numpy.sqrt(du32**2 + dv32**2) / dz32
            # get the average wind shear across the layer
            vws[i] = (vws32 + vws21) / 2.0
            continue
        elif i == 0:
            # wind shear for the top of the atmosphere (in the GRIB file)
            du = u[i,:,:] - u[i+1,:,:]
            dv = v[i,:,:] - v[i+1,:,:]
            dz = z[i,:,:] - z[i+1,:,:]
        elif i == u.shape[0]-1:
            # wind shear for the bottom of the atmosphere (in the GRIB file)
            du = u[i-1,:,:] - u[i,:,:]
            dv = v[i-1,:,:] - v[i,:,:]
            dz = z[i-1,:,:] - z[i,:,:]
        else:
            continue
        vws[i] = numpy.sqrt(du**2 + dv**2) / dz
    return vws

# constants
MPS_TO_KT = 1.94384     # meters per second to knots conversion factor
debug = False           # set this to "True" to plot out just one map

# open the grib file
rap = pygrib.open('rap/latest_rap.grib2')

# loop through the messages and get the height, u wind, and v wind messages
uwind = []
uwind_pres = []
vwind = []
vwind_pres = []
hght = []
hght_pres = []

for grb in rap:
    if grb['typeOfLevel'] == 'isobaricInhPa' and grb['name'] == 'U component of wind':
        uwind.append(grb.values)
        uwind_pres.append(grb['level'])
    elif grb['typeOfLevel'] == 'isobaricInhPa' and grb['name'] == 'V component of wind':
        vwind.append(grb.values)
        vwind_pres.append(grb['level'])
    elif grb['typeOfLevel'] == 'isobaricInhPa' and grb['name'] == 'Geopotential Height':
        hght.append(grb.values)
        hght_pres.append(grb['level'])
    else:
        continue

# get the latitude and longitude values
lats,lons = grb.latlons()

# check to make sure the pressure levels are lined up...they should be, but make sure!
assert (uwind_pres == vwind_pres == hght_pres),"Pressure levels must be the same!"

# convert everything into numpy arrays and do some more checks
uwind = numpy.array(uwind)
vwind = numpy.array(vwind)
hght = numpy.array(hght)
assert (uwind.shape == vwind.shape == hght.shape),"Wind and height fields must be the same shape!"
assert (uwind[0,:,:].shape == lats.shape),"Latitude does not match shape of vector fields!"
assert (uwind[0,:,:].shape == lons.shape),"Longitude does not match shape of vector fields!"

# run the base variables through a gaussian filter to smooth the output
sigma = 3
for i in range(uwind.shape[0]):
    uwind[i,:,:] = gfilter(uwind[i,:,:],sigma)
    vwind[i,:,:] = gfilter(vwind[i,:,:],sigma)
    hght[i,:,:] = gfilter(hght[i,:,:],sigma)

# get derived variables at each pressure level
deformation,convergence = windProducts(lons,lats,uwind,vwind)
vws = windShear(uwind,vwind,hght)
ellrod = (vws * (deformation + convergence)) * 1e7
wspd = windSpeed(uwind,vwind) * MPS_TO_KT

# create the plots of wind speed and Ellrod Index
for i in range(wspd.shape[0]):
    if debug and uwind_pres[i] != 300:
        continue
    else:
        print('%f mb pressure level' % uwind_pres[i])
        contourPlot(wspd[i,:,:],ellrod[i,:,:],uwind_pres[i])

print('Done')
