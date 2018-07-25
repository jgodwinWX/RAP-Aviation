#!/usr/bin/python3

import pygrib
import numpy
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

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
    cs_line = m.contour(x,y,var2,levels=numpy.arange(0,18000,60),colors='k',linewidths=0.5)
    cbar = m.colorbar(cs_fill,location='bottom',pad='5%')
    cbar.set_label('Wind Speed (kt)')
    plt.clabel(cs_line,cs_line.levels,inline=True,fmt='%.0f',fontsize=8)
    plt.title('RAP Analyzed %d mb Wind Speed and Geopotential' % lev)
    plt.savefig('plots/rap_wspd_%.0f.png' % lev,bbox_inches='tight')

# function for creating a basemap
def plotBasemap():
    m = Basemap(llcrnrlon=-126.138,llcrnrlat=16.281,urcrnrlon=-57.383,urcrnrlat=55.481,projection='lcc',lat_1=25.0,lat_2=25.0,lon_0=-95.0,resolution ='i')
    m.drawcoastlines(color='grey')
    m.drawstates(color='grey')
    m.drawcountries(color='grey')
    return m

# constants
MPS_TO_KT = 1.94384

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

wspd = windSpeed(uwind,vwind) * MPS_TO_KT

for i in range(wspd.shape[0]):
    contourPlot(wspd[i,:,:],hght[i,:,:],uwind_pres[i])
