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

## Set up the figure
fig = plt.figure( )
gs = gridspec.GridSpec( 2, 1, height_ratios = [ 1, 1 ] )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )

## Operations
# Connect to SQLite database.
try:
        # database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )
        # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db")
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "3.2Density_1cmSize/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

data1 = pd.read_sql( "SELECT    initial_position_x,                                         \
                                initial_position_y,                                         \
                                initial_position_z,                                         \
                                position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                velocity_x,                                                 \
                                velocity_y,                                                 \
                                velocity_z,                                                 \
                                ROUND( initial_velocity_magnitude ),                        \
                                time,                                                       \
                                crash_flag,                                                 \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( crash_flag = 1 );",                                       \
                     database )

data1.columns = [ 'init_pos_x',                                          \
                  'init_pos_y',                                          \
                  'init_pos_z',                                          \
                  'x',                                                   \
                  'y',                                                   \
                  'z',                                                   \
                  'vx',                                                  \
                  'vy',                                                  \
                  'vz',                                                  \
                  'init_vel_mag',                                        \
                  'time',                                                \
                  'crash_flag',                                          \
                  'launch_azimuth' ]

x                   = data1[ 'x' ]
y                   = data1[ 'y' ]
z                   = data1[ 'z' ]
vx                  = data1[ 'vx' ]
vy                  = data1[ 'vy' ]
vz                  = data1[ 'vz' ]
initial_velocity    = data1[ 'init_vel_mag' ]
t                   = data1[ 'time' ]
crashFlag           = data1[ 'crash_flag' ]
azimuth             = data1[ 'launch_azimuth' ]

xPositionStart      = data1[ 'init_pos_x' ]
yPositionStart      = data1[ 'init_pos_y' ]
zPositionStart      = data1[ 'init_pos_z' ]

## calculate the lat long for starting point
startRadialDistance = np.sqrt( xPositionStart**2 + yPositionStart**2 + zPositionStart**2 )
startLongitude = np.arctan2( yPositionStart, xPositionStart ) * 180.0 / np.pi
startLatitude = np.arcsin( zPositionStart / startRadialDistance ) * 180.0 / np.pi

## calculate lat long for end points
endRadialDistance = np.sqrt( x**2 + y**2 + z**2 )
endLongitude = np.arctan2( y, x ) * 180 / np.pi
endLatitude = np.arcsin( z / endRadialDistance ) * 180 / np.pi

## indicate starting point
ax1.scatter( startLongitude[0], startLatitude[0], color='black' )
ax1.text( startLongitude[0], startLatitude[0], 'start', size=20, zorder=1, color='black' )

## Plot end points on lat long map
# ax1.scatter( endLongitude, endLatitude, color='red' )
unique_reimpact_velocities = np.unique( initial_velocity )
gridsize = len( unique_reimpact_velocities )

# print gridsize
# print len(initial_velocity)
# print len(endLatitude)
# print len(endLongitude)
# sys.exit( )

hexBinPlot = ax1.hexbin( endLongitude, endLatitude,                                   \
                         C=initial_velocity,                                          \
                         cmap='jet', gridsize=gridsize )
cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax1 )

## format axis and title
ax1.set_xlabel('longitude [degree]')
ax1.set_ylabel('latitude [degree]')
cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )
ax1.set_xlim( -180.0, 180.0 )
ax1.set_yticks( np.arange( -90.0, 90.0, 10.0 ) )
ax1.set_title( 'Qualitative re-impact map' )
ax1.grid( True )

## plot the actual data points, scatter points
# ax2.scatter( endLongitude, endLatitude, s=5, c=initial_velocity, edgecolors='face', cmap='jet' )
# cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax2 )

unique_reimpact_velocities = np.unique( initial_velocity )
colors = plt.cm.Vega20( np.linspace( 0, 1, len( unique_reimpact_velocities ) ) )

for index in range( 0, len( unique_reimpact_velocities ) ):
    current_velocity = unique_reimpact_velocities[ index ]
    current_velocity_indices = np.where( initial_velocity == current_velocity )
    current_velocity_indices = current_velocity_indices[ 0 ]
    plotLatitude = endLatitude[ current_velocity_indices ]
    plotLongitude = endLongitude[ current_velocity_indices ]
    ax2.scatter( plotLongitude, plotLatitude, s=5, c=colors[ index ],                   \
                 edgecolors='face',                                                     \
                 label='$V_{launch}$ = ' + str( current_velocity ) + ' [m/s]' )

ax2.legend( ).draggable( )
ax2.set_xlabel('longitude [degree]')
ax2.set_ylabel('latitude [degree]')
# cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )
ax2.set_xlim( -180.0, 180.0 )
ax2.set_yticks( np.arange( -90.0, 90.0, 10.0 ) )
ax2.set_title( 'Actual re-impact locations' )
ax2.grid( True )

if database:
    database.close( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.suptitle( 'Regolith crash map for multiple launch velocities \n Ellipsoidal asteroid case' )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
