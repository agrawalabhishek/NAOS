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

## Operations
# Read data in csv file. data returned as a panda series.
data = pd.read_csv( '../data/ellipsoidSurfaceAcceleration.csv' )
Ux = data['Ux'].values
Uy = data['Uy'].values
Uz = data['Uz'].values
U = data['U'].values
print U
latitude = data['latitude'].values
longitude = data['longitude'].values

## Plot magnitude of acceleration due to gravity on surface of an ellipsoid using scatter map
fig = plt.figure()
ax1 = fig.add_subplot(111)
scatterPlotHeat = ax1.scatter( latitude, longitude, c=U )
cbar = plt.colorbar( scatterPlotHeat, cmap = cm.jet )
cbar.ax.set_ylabel( 'Gravitational Acceleration [m/s^2]' )
ax1.set_xlim( latitude.min(), latitude.max() )
ax1.set_ylim( longitude.min(), longitude.max() )
formatter = matplotlib.ticker.ScalarFormatter( useOffset=False )
ax1.xaxis.set_major_formatter( formatter )
ax1.yaxis.set_major_formatter( formatter )
ax1.get_yaxis().set_tick_params( direction='out' )
ax1.get_xaxis().set_tick_params( direction='out' )
ax1.set_ylabel( 'Latitude [deg]' )
ax1.set_xlabel( 'Longitude [deg]' )
ax1.set_title( 'Gravitational acceleration at ellipsoid surface (Eros)' )
# plt.ticklabel_format( style = 'sci', axis = 'x', scilimits = ( 0, 0 ) )
plt.grid()

## Plot magnitude of acceleration due to gravity on surface of an ellipsoid using contourf
fig = plt.figure()
ax1 = fig.add_subplot(111)
# find number of unique latitudes and longitudes
numberOfLatitudes = len( data['latitude'].unique() )
numberOfLongitudes = len( data['longitude'].unique() )
# make 2D arrays without changing data
y = latitude.reshape( numberOfLatitudes, numberOfLongitudes )
x = longitude.reshape( numberOfLatitudes, numberOfLongitudes )
z = U.reshape( numberOfLatitudes, numberOfLongitudes )
contourHeatPlot = plt.contourf( y, x, z, cmap=cm.jet )
cbar = plt.colorbar( contourHeatPlot, cmap=cm.jet )
cbar.ax.set_ylabel( 'Gravitational Acceleration [m/s^2]' )
ax1.set_xlim( latitude.min(), latitude.max() )
ax1.set_ylim( longitude.min(), longitude.max() )
formatter = matplotlib.ticker.ScalarFormatter( useOffset=False )
ax1.xaxis.set_major_formatter( formatter )
ax1.yaxis.set_major_formatter( formatter )
ax1.get_yaxis().set_tick_params( direction='out' )
ax1.get_xaxis().set_tick_params( direction='out' )
ax1.set_ylabel( 'Latitude [deg]' )
ax1.set_xlabel( 'Longitude [deg]' )
ax1.set_title( 'Gravitational acceleration at ellipsoid surface (Eros)' )
# plt.ticklabel_format( style = 'sci', axis = 'x', scilimits = ( 0, 0 ) )
plt.grid()
plt.show()

# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# heatmap = plt.pcolor( x, y, z, cmap=cm.jet )
# plt.show()

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
