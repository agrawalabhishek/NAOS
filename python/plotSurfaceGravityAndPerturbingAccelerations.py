'''
Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
Distributed under the MIT License.
See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
'''

# Set up modules and packages
# I/O
import csv
from pprint import pprint

# Numerical
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import math

# System
import sys
import time
from tqdm import tqdm

# Get plotting packages
import matplotlib
import matplotlib.colors
import matplotlib.axes
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import rcParams
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.tri as tri

print ""
print "---------------------------------------------------------------------------------"
print "                                 NAOS                                            "
print "                                                                                 "
print "         Copyright (c) 2016, A. Agrawal (abhishek.agrawal@protonmail.com)        "
print "---------------------------------------------------------------------------------"
print ""

# Start timer.
start_time = time.time( )

def fmt(x, pos):
    a, b = '{:.7e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

## Operations
# Read data in csv file. data returned as a panda series.
data = pd.read_csv( '../data/surface_gravity_and_perturbing_acceleration/equatorial_and_circular_case/allAccelerations.csv' )
Ux = data['xGravAcceleration'].values
Uy = data['yGravAcceleration'].values
Uz = data['zGravAcceleration'].values
U = data['gravAcceleration'].values

sunPeriapsisTBP = data['thirdBodyAcceleration_sunPeriapsis'].values
sunPeriapsisSRP = data['srpAcceleration_sunPeriapsis'].values

sunApoapsisTBP = data['thirdBodyAcceleration_sunApoapsis'].values
sunApoapsisSRP = data['srpAcceleration_sunApoapsis'].values

latitude = data['latitude'].values
longitude = data['longitude'].values

# find number of unique latitudes and longitudes
numberOfLatitudes = len( data['latitude'].unique() )
numberOfLongitudes = len( data['longitude'].unique() )

## Plot magnitude of acceleration due to gravity on surface of an ellipsoid using contourf
fig = plt.figure()
ax1 = fig.add_subplot(111)
# make 2D arrays without changing data
y = latitude.reshape( numberOfLatitudes, numberOfLongitudes )
x = longitude.reshape( numberOfLatitudes, numberOfLongitudes )
z = U.reshape( numberOfLatitudes, numberOfLongitudes )
contourHeatPlot = plt.contourf( x, y, z, cmap=cm.jet )
cbar = plt.colorbar( contourHeatPlot, cmap=cm.jet )
cbar.ax.set_ylabel( 'Gravitational Acceleration [m/s^2]' )
ax1.set_ylim( latitude.min(), latitude.max() )
ax1.set_xlim( longitude.min(), longitude.max() )
formatter = matplotlib.ticker.ScalarFormatter( useOffset=False )
ax1.xaxis.set_major_formatter( formatter )
ax1.yaxis.set_major_formatter( formatter )
ax1.get_yaxis().set_tick_params( direction='out' )
ax1.get_xaxis().set_tick_params( direction='out' )
ax1.set_ylabel( 'Latitude [deg]' )
ax1.set_xlabel( 'Longitude [deg]' )
ax1.set_title( 'Gravitational acceleration at ellipsoid surface (Eros)' )

plt.grid()
plt.show()

## Plot magnitude of acceleration due to sun's third body effect on surface of an ellipsoid using contourf
fig = plt.figure()
ax1 = fig.add_subplot(111)
# make 2D arrays without changing data
y = latitude.reshape( numberOfLatitudes, numberOfLongitudes )
x = longitude.reshape( numberOfLatitudes, numberOfLongitudes )
z = sunPeriapsisTBP.reshape( numberOfLatitudes, numberOfLongitudes )
contourHeatPlot = plt.contourf( x, y, z, cmap=cm.jet )
cbar = plt.colorbar( contourHeatPlot, cmap=cm.jet )
cbar.ax.set_ylabel( 'Sun Third Body Perturbing Acceleration [m/s^2]' )
ax1.set_ylim( latitude.min(), latitude.max() )
ax1.set_xlim( longitude.min(), longitude.max() )
formatter = matplotlib.ticker.ScalarFormatter( useOffset=False )
ax1.xaxis.set_major_formatter( formatter )
ax1.yaxis.set_major_formatter( formatter )
ax1.get_yaxis().set_tick_params( direction='out' )
ax1.get_xaxis().set_tick_params( direction='out' )
ax1.set_ylabel( 'Latitude [deg]' )
ax1.set_xlabel( 'Longitude [deg]' )
# ax1.set_title( 'Third body effect at Eros surface (Sun Periapsis)' )
ax1.set_title( 'Third body effect at Eros surface' )

## Plot magnitude of acceleration due to srp effect on surface of an ellipsoid using contourf
fig = plt.figure()
ax1 = fig.add_subplot(111)
# make 2D arrays without changing data
y = latitude.reshape( numberOfLatitudes, numberOfLongitudes )
x = longitude.reshape( numberOfLatitudes, numberOfLongitudes )
z = sunPeriapsisSRP.reshape( numberOfLatitudes, numberOfLongitudes )
contourHeatPlot = plt.contourf( x, y, z, cmap=cm.jet )
cbar = plt.colorbar( contourHeatPlot, cmap=cm.jet, format=matplotlib.ticker.FuncFormatter(fmt) )
cbar.ax.set_ylabel( 'Solar Radiation Pressure Acceleration [m/s^2]' )
ax1.set_ylim( latitude.min(), latitude.max() )
ax1.set_xlim( longitude.min(), longitude.max() )
# formatter = matplotlib.ticker.ScalarFormatter( useOffset=False )
ax1.xaxis.set_major_formatter( formatter )
ax1.yaxis.set_major_formatter( formatter )
ax1.get_yaxis().set_tick_params( direction='out' )
ax1.get_xaxis().set_tick_params( direction='out' )
ax1.set_ylabel( 'Latitude [deg]' )
ax1.set_xlabel( 'Longitude [deg]' )
# ax1.set_title( 'SRP acceleration at Eros surface (Sun Periapsis)' )
ax1.set_title( 'SRP acceleration at Eros surface' )

## Plot magnitude of acceleration due to sun's third body effect on surface of an ellipsoid using contourf
fig = plt.figure()
ax1 = fig.add_subplot(111)
# make 2D arrays without changing data
y = latitude.reshape( numberOfLatitudes, numberOfLongitudes )
x = longitude.reshape( numberOfLatitudes, numberOfLongitudes )
z = sunApoapsisTBP.reshape( numberOfLatitudes, numberOfLongitudes )
contourHeatPlot = plt.contourf( x, y, z, cmap=cm.jet )
cbar = plt.colorbar( contourHeatPlot, cmap=cm.jet )
cbar.ax.set_ylabel( 'Sun Third Body Perturbing Acceleration [m/s^2]' )
ax1.set_ylim( latitude.min(), latitude.max() )
ax1.set_xlim( longitude.min(), longitude.max() )
formatter = matplotlib.ticker.ScalarFormatter( useOffset=False )
ax1.xaxis.set_major_formatter( formatter )
ax1.yaxis.set_major_formatter( formatter )
ax1.get_yaxis().set_tick_params( direction='out' )
ax1.get_xaxis().set_tick_params( direction='out' )
ax1.set_ylabel( 'Latitude [deg]' )
ax1.set_xlabel( 'Longitude [deg]' )
ax1.set_title( 'Third body effect at Eros surface (Sun Apoapsis)' )

## Plot magnitude of acceleration due to srp effect on surface of an ellipsoid using contourf
fig = plt.figure()
ax1 = fig.add_subplot(111)
# make 2D arrays without changing data
y = latitude.reshape( numberOfLatitudes, numberOfLongitudes )
x = longitude.reshape( numberOfLatitudes, numberOfLongitudes )
z = sunApoapsisSRP.reshape( numberOfLatitudes, numberOfLongitudes )
contourHeatPlot = plt.contourf( x, y, z, cmap=cm.jet )
cbar = plt.colorbar( contourHeatPlot, cmap=cm.jet, format=matplotlib.ticker.FuncFormatter(fmt) )
cbar.ax.set_ylabel( 'Solar Radiation Pressure Acceleration [m/s^2]' )
ax1.set_ylim( latitude.min(), latitude.max() )
ax1.set_xlim( longitude.min(), longitude.max() )
# formatter = matplotlib.ticker.ScalarFormatter( useOffset=False )
ax1.xaxis.set_major_formatter( formatter )
ax1.yaxis.set_major_formatter( formatter )
ax1.get_yaxis().set_tick_params( direction='out' )
ax1.get_xaxis().set_tick_params( direction='out' )
ax1.set_ylabel( 'Latitude [deg]' )
ax1.set_xlabel( 'Longitude [deg]' )
ax1.set_title( 'SRP acceleration at Eros surface (Sun Apoapsis)' )

plt.grid()
plt.show()

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
