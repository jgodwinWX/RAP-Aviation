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
def contourPlot(var1,var2,lev,title,save):
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
    if save == 'ellrod':
        values = numpy.arange(4,40,2)
        fmat = '%.0f'
    elif save == 'ri':
        values = numpy.arange(0,1,0.1)
        fmat = '%.1f'
    else:
        values = numpy.arange(numpy.nanmin(var2),numpy.nanmax(var2),(numpy.nanmax(var2)-numpy.nanmin(var2))/10.0)
        fmat = '%f'
    cs_line = m.contour(x,y,var2,levels=values,colors='r',linewidths=0.5)
    cbar = m.colorbar(cs_fill,location='bottom',pad='5%')
    cbar.set_label('Wind Speed (kt)')
    plt.clabel(cs_line,cs_line.levels,inline=True,fmt=fmat,fontsize=6)
    plt.title('RAP Analyzed %d mb %s' % (lev,title))
    plt.savefig('plots/rap_%s_%.0f.png' % (save,lev),bbox_inches='tight')

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

def richardsonNumber(theta,z,vws):
    GRAVITY = 9.80665   # acceleration of gravity (m**2/s**2)
    # first we compute the theta lapse rate
    ri = numpy.zeros(theta.shape)
    for i in range(theta.shape[0]):
        if i > 0 and i < theta.shape[0] - 1:
            dt21 = theta[i-1,:,:] - theta[i,:,:]
            dz21 = z[i-1,:,:] - z[i,:,:]
            theta_lr21 = dt21 / dz21
            dt32 = theta[i,:,:] - theta[i+1,:,:]
            dz32 = z[i,:,:] - z[i+1,:,:]
            theta_lr32 = dt32 / dz32
            theta_lr = (theta_lr32 + theta_lr21) / 2.0
        elif i == 0:
            dt = theta[i,:,:] - theta[i+1,:,:]
            dz = z[i,:,:] - z[i+1,:,:]
            theta_lr = dt / dz
        elif i == theta.shape[0]-1:
            dt = theta[i-1,:,:] - theta[i,:,:]
            dz = z[i-1,:,:] - z[i,:,:]
            theta_lr = dt / dz
        else:
            continue
        # Brunt-Vaisala Frequency
        bvf = numpy.sqrt((GRAVITY / theta[i,:,:]) * (theta_lr))
        
        # actual Richardson Number computation
        ri[i] = bvf**2 / vws[i,:,:]**2
    return ri

# function for computing potential temperature
def potentialTemperature(t,p):
    KAPPA = 0.2854  # Poisson Constant
    theta = numpy.zeros(t.shape)
    for i in range(t.shape[0]):
        theta[i,:,:] = t[i,:,:] * (1000.0 / p[i])**KAPPA
    return theta

# constants
MPS_TO_KT = 1.94384     # meters per second to knots conversion factor
sigma = 3               # smoothing factor used by Gaussian Filter
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
temp = []
temp_pres = []

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
    elif grb['typeOfLevel'] == 'isobaricInhPa' and grb['name'] == 'Temperature':
        temp.append(grb.values)
        temp_pres.append(grb['level'])
    else:
        continue

# get the latitude and longitude values
lats,lons = grb.latlons()

# check to make sure the pressure levels are lined up...they should be, but make sure!
assert (uwind_pres == vwind_pres == hght_pres == temp_pres),"Pressure levels must be the same!"

# convert everything into numpy arrays and do some more checks
uwind = numpy.array(uwind)
vwind = numpy.array(vwind)
hght = numpy.array(hght)
temp = numpy.array(temp)
assert (uwind.shape == vwind.shape == hght.shape == temp.shape),"All fields must be the same shape!"
assert (uwind[0,:,:].shape == lats.shape),"Latitude does not match shape of vector fields!"
assert (uwind[0,:,:].shape == lons.shape),"Longitude does not match shape of vector fields!"

# run the base variables through a gaussian filter to smooth the output
for i in range(uwind.shape[0]):
    uwind[i,:,:] = gfilter(uwind[i,:,:],sigma)
    vwind[i,:,:] = gfilter(vwind[i,:,:],sigma)
    hght[i,:,:] = gfilter(hght[i,:,:],sigma)
    temp[i,:,:] = gfilter(temp[i,:,:],sigma)

# get derived variables at each pressure level
deformation,convergence = windProducts(lons,lats,uwind,vwind)
vws = windShear(uwind,vwind,hght)
ellrod = (vws * (deformation + convergence)) * 1e7
wspd = windSpeed(uwind,vwind) * MPS_TO_KT
theta = potentialTemperature(temp,temp_pres)
ri = richardsonNumber(theta,hght,vws)

# create the plots of wind speed and Ellrod Index
for i in range(wspd.shape[0]):
    if debug and uwind_pres[i] != 300:
        continue
    else:
        print('%.0f mb pressure level' % uwind_pres[i])
        contourPlot(wspd[i,:,:],ellrod[i,:,:],uwind_pres[i],'Wind Speed and Ellrod Index','ellrod')
        contourPlot(wspd[i,:,:],ri[i,:,:],uwind_pres[i],'Wind Speed and Richardson Number','ri')

print('Done')
