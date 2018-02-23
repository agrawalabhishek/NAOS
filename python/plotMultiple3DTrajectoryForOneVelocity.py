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
## data in body frame
# data = pd.read_csv( '../data/trajectory_for_different_launch_azimuth/regolithTrajectoryAtAzimuth210.csv' )
# x = data[ 'x' ].values
# y = data[ 'y' ].values
# z = data[ 'z' ].values
# vx = data[ 'vx' ].values
# vy = data[ 'vy' ].values
# vz = data[ 'vz' ].values
# t = data[ 't' ].values

## Plot the ellipsoidal shape of the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593

# Connect to SQLite database.
try:
        # database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

azimuthRange = tuple( ( np.arange( 0.0, 360.0, 20.0 ) ).tolist( ) )

data1 = pd.read_sql( "SELECT     trajectory_id                                                   \
                      FROM       regolith_trajectory_results                                     \
                      WHERE      ROUND( launch_azimuth ) IN " + str(azimuthRange) + "            \
                      AND        ROUND( initial_velocity_magnitude ) = 5.0                       \
                      AND        crash_flag = 1;",                                               \
                      database )

data1.columns = [ 'trajectory_id' ]
trajectory_id = data1[ 'trajectory_id' ]
trajectory_id = tuple( trajectory_id.tolist( ) )

data = pd.read_sql( "SELECT     position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                ROUND( initial_velocity_magnitude ),                            \
                                inertial_position_x,                                            \
                                inertial_position_y,                                            \
                                inertial_position_z,                                            \
                                ROUND( launch_azimuth ),                                        \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      trajectory_id IN " + str( trajectory_id ) + ";",                \
                     database )

if database:
    database.close( )

data.columns = [ 'x',                                                                           \
                 'y',                                                                           \
                 'z',                                                                           \
                 'velocity_magnitude',                                                          \
                 'inertial_x',                                                                  \
                 'inertial_y',                                                                  \
                 'inertial_z',                                                                  \
                 'launch_azimuth',                                                              \
                 'time' ]

x                   = data[ 'x' ]
y                   = data[ 'y' ]
z                   = data[ 'z' ]
velocityMagnitude   = data[ 'velocity_magnitude' ]
inertial_x          = data[ 'inertial_x' ]
inertial_y          = data[ 'inertial_y' ]
inertial_z          = data[ 'inertial_z' ]
launchAzimuth       = data[ 'launch_azimuth' ]
t                   = data[ 'time' ]

print "Processing data now..."

## Set up the figure
fig = plt.figure( )
gs = gridspec.GridSpec( 2, 1, height_ratios = [ 1, 1 ] )
ax1 = plt.subplot( gs[ 0 ], projection = '3d' )
ax2 = plt.subplot( gs[ 1 ], projection = '3d' )

plt.suptitle( 'Particle trajectory around asteroid Eros (Body frame) \n'
               + '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '
               + 'Longest edge, Reimpact case' )

## draw the ellipsoid
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

newColor = colors.cnames["slategray"]
surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         rstride=5, cstride=5, alpha=0.5 )
surf.set_facecolor( newColor )
surf.set_linewidth( 0.1 )

surf = ax2.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         rstride=5, cstride=5, alpha=0.5 )
surf.set_facecolor( newColor )
surf.set_linewidth( 0.1 )

## unique azimuth
unique_azimuths = np.unique( launchAzimuth )
colors = plt.cm.Vega20( np.linspace( 0, 1, len( unique_azimuths ) ) )

## Plot 3D trajectory of the orbiting particle
for index in range( 0, len( unique_azimuths ) ):
    current_azimuth = unique_azimuths[ index ]
    data_indices = np.where( launchAzimuth == current_azimuth )
    data_indices = data_indices[ 0 ]
    xPlot = (x[ data_indices ]).tolist( )
    yPlot = (y[ data_indices ]).tolist( )
    zPlot = (z[ data_indices ]).tolist( )
    ax1.plot( xPlot, yPlot, zPlot, zdir = 'z', c=colors[index],                \
              label='Launch azimuth = ' + str( current_azimuth ) + '[deg]' )

ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax1.set_title( 'Body Frame' )
ax1.grid( True )
# ax1.legend( ).draggable( )
plt.show( )
sys.exit( )

## plot the 3d trajectory in inertial frame
ax2.plot( inertial_x, inertial_y, inertial_z, zdir = 'z', color=colors.cnames["purple"] )

ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.set_zlabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax2.set_title( 'Inertial Frame' )
ax2.grid( True )

if database:
    database.close( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## show the plot
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
