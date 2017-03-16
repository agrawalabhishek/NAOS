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
gs = gridspec.GridSpec( 1, 3 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

## Operations
# Connect to SQLite database.
try:
        # database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db")

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

# start position for all cases is the same
xPositionStart = xPositionStart[ 0 ]
yPositionStart = yPositionStart[ 0 ]
zPositionStart = zPositionStart[ 0 ]

unique_reimpact_velocities = np.unique( initial_velocity )
colors = plt.cm.Vega20( np.linspace( 0, 1, len( unique_reimpact_velocities ) ) )

for index in range( 0, len( unique_reimpact_velocities ) ):
    current_velocity = unique_reimpact_velocities[ index ]
    current_velocity_indices = np.where( initial_velocity == current_velocity )
    current_velocity_indices = current_velocity_indices[ 0 ]

    azimuthPlot = azimuth[ current_velocity_indices ]
    xPlot = x[ current_velocity_indices ]
    yPlot = y[ current_velocity_indices ]
    zPlot = z[ current_velocity_indices ]
    reimpactRange = np.sqrt( (xPlot - xPositionStart)**2 + (yPlot - yPositionStart)**2 + (zPlot - zPositionStart)**2 )

    if current_velocity <= 6.0:
        ax1.scatter( azimuthPlot, reimpactRange, s=5, c=colors[ index ],                   \
                     edgecolors='face',                                                    \
                     label='$V_{launch}$ = ' + str( current_velocity ) + ' [m/s]' )
    elif current_velocity > 6.0 and current_velocity <= 11.0:
        ax2.scatter( azimuthPlot, reimpactRange, s=5, c=colors[ index ],                   \
                     edgecolors='face',                                                    \
                     label='$V_{launch}$ = ' + str( current_velocity ) + ' [m/s]' )
    elif current_velocity > 11.0:
        ax3.scatter( azimuthPlot, reimpactRange, s=5, c=colors[ index ],                   \
                     edgecolors='face',                                                    \
                     label='$V_{launch}$ = ' + str( current_velocity ) + ' [m/s]' )

ax1.legend( ).draggable( )
ax1.set_xlabel('Launch azimuth [deg]')
ax1.set_ylabel('Range [m]')
ax1.grid( True )

ax2.legend( ).draggable( )
ax2.set_xlabel('Launch azimuth [deg]')
ax2.set_ylabel('Range [m]')
ax2.grid( True )

ax3.legend( ).draggable( )
ax3.set_xlabel('Launch azimuth [deg]')
ax3.set_ylabel('Range [m]')
ax3.grid( True )

if database:
    database.close( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.suptitle( 'Reimpact point range from launch location' )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
