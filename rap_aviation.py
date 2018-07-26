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
    elif save == 'dt':
        values = numpy.arange(5,100,5)
        fmat = '%.0f'
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

def importVars(grib):
    uwind = []
    uwind_pres = []
    vwind = []
    vwind_pres = []
    hght = []
    hght_pres = []
    temp = []
    temp_pres = []

    for grb in grib:
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
    lats,lons = grb.latlons()
    return uwind,uwind_pres,vwind,vwind_pres,hght,hght_pres,temp,temp_pres,lats,lons

# check the pressure levels
def pressureCheck(var1,var2,var3,var4):
    assert (var1 == var2 == var3 == var4),"Pressure levels must be the same!"

# check the shapes
def shapeCheck(var1,var2,var3,var4,la,lo):
    assert (var1.shape == var2.shape == var3.shape == var4.shape),"All fields must be the same shape!"
    assert (var1[0,:,:].shape == la.shape),"Latitude does not match shape of vector fields!"
    assert (var1[0,:,:].shape == lo.shape),"Longitude does not match shape of vector fields!"

# constants
MPS_TO_KT = 1.94384     # meters per second to knots conversion factor
sigma = 3               # smoothing factor used by Gaussian Filter
debug = False           # set this to "True" to plot out just one map

# open the grib file
rap = pygrib.open('rap/latest_rap.grib2')
fcst = pygrib.open('rap/rap_fcst03.grib2')

uwind_now,uwind_pres_now,vwind_now,vwind_pres_now,hght_now,hght_pres_now,temp_now,temp_pres_now,lats,lons = importVars(rap)
uwind_fcst,uwind_pres_fcst,vwind_fcst,vwind_pres_fcst,hght_fcst,hght_pres_fcst,temp_fcst,temp_pres_fcst,lats_fcst,lons_fcst = importVars(fcst)

# make sure the pressure levels line up
pressureCheck(uwind_pres_now,vwind_pres_now,hght_pres_now,temp_pres_now)
pressureCheck(uwind_pres_fcst,vwind_pres_fcst,hght_pres_fcst,temp_pres_fcst)

# convert everything into numpy arrays and do some more checks
uwind_now = numpy.array(uwind_now)
vwind_now = numpy.array(vwind_now)
hght_now = numpy.array(hght_now)
temp_now = numpy.array(temp_now)
uwind_fcst = numpy.array(uwind_fcst)
vwind_fcst = numpy.array(vwind_fcst)
hght_fcst = numpy.array(hght_fcst)
temp_fcst = numpy.array(temp_fcst)

# check the shapes of everything
shapeCheck(uwind_now,vwind_now,hght_now,temp_now,lats,lons)
shapeCheck(uwind_fcst,vwind_fcst,hght_fcst,temp_fcst,lats,lons)

# run the base variables through a gaussian filter to smooth the output
for i in range(uwind_now.shape[0]):
    uwind_now[i,:,:] = gfilter(uwind_now[i,:,:],sigma)
    vwind_now[i,:,:] = gfilter(vwind_now[i,:,:],sigma)
    hght_now[i,:,:] = gfilter(hght_now[i,:,:],sigma)
    temp_now[i,:,:] = gfilter(temp_now[i,:,:],sigma)
    uwind_fcst[i,:,:] = gfilter(uwind_fcst[i,:,:],sigma)
    vwind_fcst[i,:,:] = gfilter(vwind_fcst[i,:,:],sigma)
    hght_fcst[i,:,:] = gfilter(hght_fcst[i,:,:],sigma)
    temp_fcst[i,:,:] = gfilter(temp_fcst[i,:,:],sigma)

# get derived variables at each pressure level
deformation_now,convergence_now = windProducts(lons,lats,uwind_now,vwind_now)
deformation_fcst,convergence_fcst = windProducts(lons,lats,uwind_fcst,vwind_fcst)
vws_now = windShear(uwind_now,vwind_now,hght_now)
ellrod_now = (vws_now * (deformation_now + convergence_now)) * 1e7
wspd_now = windSpeed(uwind_now,vwind_now) * MPS_TO_KT
theta_now = potentialTemperature(temp_now,temp_pres_now)
ri_now = richardsonNumber(theta_now,hght_now,vws_now)
div_tend = (-convergence_fcst + convergence_now) / (3600*3)

# create the plots of wind speed and Ellrod Index
for i in range(wspd_now.shape[0]):
    if debug and uwind_pres_now[i] != 300:
        continue
    else:
        print('%.0f mb pressure level' % uwind_pres_now[i])
        contourPlot(wspd_now[i,:,:],ellrod_now[i,:,:],uwind_pres_now[i],'Wind Speed and Ellrod Index','ellrod')
        contourPlot(wspd_now[i,:,:],ri_now[i,:,:],uwind_pres_now[i],'Wind Speed and Richardson Number','ri')
        contourPlot(wspd_now[i,:,:],div_tend[i,:,:]*1e9,uwind_pres_now[i],'Wind Speed and Divergence Tendency','dt')

print('Done')
