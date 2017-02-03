'''
Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
Distributed under the MIT License.
See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
'''

# Set up modules and packages
# I/O
import csv
from pprint import pprint
import sqlite3

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
import matplotlib.colors as colors
import matplotlib.axes
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import rcParams
from matplotlib import cm
from matplotlib import gridspec
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
data = pd.read_csv( '../data/sun_asteroid_2BP/equatorial_and_circular_case/sunAsteroid2BP.csv' )

xBody           = data["x_body_frame"].values
yBody           = data["y_body_frame"].values
zBody           = data["z_body_frame"].values
vxBody          = data["vx_body_frame"].values
vyBody          = data["vy_body_frame"].values
vzBody          = data["vz_body_frame"].values
t               = data["t"].values
xInertial       = data["x_inertial_frame"].values
yInertial       = data["y_inertial_frame"].values
zInertial       = data["z_inertial_frame"].values
vxInertial      = data["vx_inertial_frame"].values
vyInertial      = data["vy_inertial_frame"].values
vzInertial      = data["vz_inertial_frame"].values
sma             = data["sma"].values
eccentricity    = data["eccentricity"].values
inclination     = data["inclination"].values
raan            = data["raan"].values
aop             = data["aop"].values
ta              = data["ta"].values
eccentricAnomaly= data["eccentric_anomaly"].values
meanAnomaly     = data["mean_anomaly"].values
jacobian        = data["jacobi"].values
energy          = data["energy"].values

radiusPerigee = sma * ( 1.0 - eccentricity )
print "radius Perigee = " + str( radiusPerigee[ 0 ] )

radiusApogee = sma * ( 1.0 + eccentricity )
print "radius Apogee = " + str( radiusApogee[ 0 ] )

# extract part of the data since the sun ephemeris comprises of multiple revolutions
# in days to slice, the divide by ten option is because the data is only saved every 10 seconds in
# the input csv file. so index value is = desired_time / 10
# timeToSlice = 150 * 24 * 60 * 60 / 10
# xInertialNew = xInertial[ 0:timeToSlice ]
# yInertialNew = yInertial[ 0:timeToSlice ]
# zInertialNew = zInertial[ 0:timeToSlice ]
# tNew         = t[ 0:timeToSlice ]

sunRadialDistance = np.sqrt( xInertial**2 + yInertial**2 + zInertial**2 )
print "min. radial distance of sun from the asteroid =  " + str( min( sunRadialDistance ) )
print "max. radial distance of sun from the asteroid =  " + str( max( sunRadialDistance ) )

oneAstronomicalUnit = 149597870700.0
fig = plt.figure( )
ax1 = fig.add_subplot( 111 )
ax1.plot( t/(24.0*60.0*60.0), ( sunRadialDistance/oneAstronomicalUnit ) )
# ax1.plot( t/(24.0*60.0*60.0), radiusPerigee/oneAstronomicalUnit, color=colors.cnames['purple'], label='Periapsis distance' )
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
ax1.set_xlabel( "time [Earth days]" )
ax1.set_ylabel( "Radial distance to Sun [AU]" )
# plt.legend( )
plt.grid( )
plt.show( )

fig = plt.figure( )
ax1 = fig.add_subplot( 311 )
ax1.plot( t/(24.0*60.0*60.0), ta )
ax1.set_xlabel( "time [Earth days]" )
ax1.set_ylabel( "True anomaly [deg]" )
ax1.grid( )

ax2 = fig.add_subplot( 312 )
ax2.plot( t/(24.0*60.0*60.0), eccentricAnomaly%360.0 )
ax2.set_xlabel( "time [Earth days]" )
ax2.set_ylabel( "Eccentric anomaly [deg]" )
ax2.grid( )

ax3 = fig.add_subplot( 313 )
ax3.plot( t/(24.0*60.0*60.0), meanAnomaly%360.0 )
ax3.set_xlabel( "time [Earth days]" )
ax3.set_ylabel( "Mean anomaly [deg]" )
ax3.grid( )

plt.show( )

## find array index for farthest approach of sun to the asteroid
# farthestApproachIndex = np.where( sunRadialDistance == max( sunRadialDistance ) )
# # farthestApproachIndex = np.where( sunRadialDistance == radiusApogee[ 0 ] )
# print "farthest approach index = " + str( farthestApproachIndex )
# timeOfFarthestApproach = t[ farthestApproachIndex ]

# print "time of farthest approach = " + str( timeOfFarthestApproach )

## error in solving the kepler problem
apogeeError = np.abs( max( sunRadialDistance ) - radiusApogee[ 0 ] )
print "error in apogee calculation = " + str( apogeeError / 1000.0 )

perigeeError = np.abs( min( sunRadialDistance ) - radiusPerigee[ 0 ] )
print "error in perigee calculation = " + str( perigeeError / 1000.0 )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
