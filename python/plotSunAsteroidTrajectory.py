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
data = pd.read_csv( '../data/sun_asteroid_2BP/equatorial_elliptical_TA0deg/sunAsteroid2BP.csv' )

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
jacobian        = data["jacobi"].values
energy          = data["energy"].values

# convert time in seconds to earth days
t = t / ( 24.0 * 60.0 * 60.0 )

## Set up the figure
fig = plt.figure( )
ax1 = fig.add_subplot( 111, projection = '3d' )

## plot the Trajectory
ax1.plot( xInertial, yInertial, zInertial, color=colors.cnames['purple'] )

ax1.scatter( xInertial[0], yInertial[0], zInertial[0], s=200,                                  \
             c=colors.cnames['orange'], marker=(5,1), label='Sun' )

ax1.scatter( 0.0, 0.0, 0.0, s=100, c=colors.cnames['black'], marker=(5,0), label='Eros' )

ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')
ax1.set_title('Asteroid-centric orbit trajectory of Sun')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax1.legend( ).draggable( )
ax1.grid(True)

## Plot trajectory in xy plane (inertial frame)
fig = plt.figure( )
ax1 = fig.add_subplot( 111 )

ax1.plot( xInertial, yInertial, color=colors.cnames['purple'] )

ax1.scatter( xInertial[ 0 ], yInertial[ 0 ], s=200,                                            \
             c=colors.cnames['orange'], marker=(5,1), label='Sun' )

ax1.scatter( 0.0, 0.0, s=100, c=colors.cnames['black'], marker=(5,0), label='Eros' )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_title('Sun trajectory around the asteroid (Inertial frame)')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax1.legend( ).draggable( )
ax1.grid(True)
ax1.axis('equal')

## Plot trajectory in xy plane (body frame)
# fig = plt.figure( )
# ax1 = fig.add_subplot( 111 )

# numberOfDaysToPlot = 0.5
# secondsToPlot = numberOfDaysToPlot * 24.0 * 60.0 * 60.0
# dataSaveInterval = 50.0
# finalIndex = secondsToPlot / dataSaveInterval
# indexRange = range( 0, int( finalIndex ) )

# ax1.plot( xBody[ indexRange ], yBody[ indexRange ], color=colors.cnames['purple'] )

# ax1.scatter( xBody[ 0 ], yBody[ 0 ], s=200,                                            \
#              c=colors.cnames['orange'], marker=(5,1), label='Sun' )

# ax1.scatter( 0.0, 0.0, s=100, c=colors.cnames['black'], marker=(5,0), label='Eros' )
# ax1.set_xlabel('x [m]')
# ax1.set_ylabel('y [m]')
# ax1.set_title('Sun trajectory around the asteroid (Body frame)')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
# ax1.legend( ).draggable( )
# ax1.grid(True)
# ax1.axis('equal')

## show plot
plt.show( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
