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
alpha = 20000.0
beta = 7000.0
gamma = 7000.0

## Operations
# Connect to SQLite database.
try:
    database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                               + "multiple_launch_velocity/"
                               + "simulation_time_9_months/"
                               + "longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

azimuthRange = np.arange( 0.0, 360.0, 5.0 ).tolist( )
azimuthRange = tuple( azimuthRange )

data1 = pd.read_sql( "SELECT    initial_position_x,                                             \
                                initial_position_y,                                             \
                                initial_position_z,                                             \
                                position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth )                                         \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ( crash_flag = 1 )                                              \
                     AND        ROUND( initial_velocity_magnitude ) = 4.0                       \
                     AND        ROUND( launch_azimuth ) IN " + str( azimuthRange ) + ";",       \
                     database )

data1.columns = [ 'initial_position_x',                                     \
                  'initial_position_y',                                     \
                  'initial_position_z',                                     \
                  'x',                                                      \
                  'y',                                                      \
                  'z',                                                      \
                  'init_vel_mag',                                           \
                  'launch_azimuth' ]

x                   = data1[ 'x' ]
y                   = data1[ 'y' ]
z                   = data1[ 'z' ]
initial_velocity    = data1[ 'init_vel_mag' ]
azimuth             = data1[ 'launch_azimuth' ]

xPositionStart      = data1[ 'initial_position_x' ]
yPositionStart      = data1[ 'initial_position_y' ]
zPositionStart      = data1[ 'initial_position_z' ]

## calculate the lat long for starting point
startRadialDistance = np.sqrt( xPositionStart**2 + yPositionStart**2 + zPositionStart**2 )
startLongitude = np.arctan2( yPositionStart, xPositionStart ) * 180.0 / np.pi
startLatitude = np.arcsin( zPositionStart / startRadialDistance ) * 180.0 / np.pi

## calculate lat long for end points
endRadialDistance = np.sqrt( x**2 + y**2 + z**2 )
endLongitude = np.arctan2( y, x ) * 180 / np.pi
endLatitude = np.arcsin( z / endRadialDistance ) * 180 / np.pi

## calculate range to reimpact point from launch point
xComp = xPositionStart - x
yComp = yPositionStart - y
zComp = zPositionStart - z
crashRange = np.sqrt( xComp**2 + yComp**2 + zComp**2 )

ax1.plot( azimuth, crashRange,
          marker='o',
          markersize=4.0,
          linewidth=0.5,
          linestyle='dotted' )

ax1.grid(True)
ax1.set_xlabel( 'Launch azimuth [deg]' )
ax1.set_ylabel( 'Range [m]' )
ax1.set_xticks( np.arange( 0.0, 355.0, 10.0 ) )
plt.xticks( rotation='vertical' )

if database:
    database.close( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.suptitle( 'Range to reimpact point from launch location for V$_{launch}$ = ' + str( initial_velocity[ 0 ] )
              + '[m/s] \n Ellipsoid longest edge' )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
