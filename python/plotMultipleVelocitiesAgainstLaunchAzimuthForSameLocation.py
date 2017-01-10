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
Wz = 0.00033118202125129593
mu = 876514

maxTime = 1.0 * 30.0 * 24.0 * 60.0 * 60.0

## Set up the figure
fig = plt.figure( )
gs = gridspec.GridSpec( 2, 1, height_ratios = [ 1, 1 ] )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
# ax1 = fig.add_subplot(111)
# ax2 = fig.add_subplot(212)

## Operations
# Connect to SQLite database.
try:
        database = sqlite3.connect("../data/regolith_launched_from_leading_edge/leadingEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

## get the local directional escape speeds in rotating and inertial frame for each launch azimuth
data4 = pd.read_sql( "SELECT    directional_escape_speed,                                   \
                                directional_inertial_escape_speed,                          \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( start_flag = 1 );",                                       \
                     database )

data4.columns = [ 'directional_escape_speed',                                               \
                  'directional_inertial_escape_speed',                                      \
                  'directional_escape_azimuth']

directionalEscapeSpeed          = data4[ 'directional_escape_speed' ]
inertialDirectionalEscapeSpeed  = data4[ 'directional_inertial_escape_speed' ]
directionalEscapeAzimuth        = data4[ 'directional_escape_azimuth' ]

## get data for escape cases
data1 = pd.read_sql( "SELECT    initial_velocity_x,                                         \
                                initial_velocity_y,                                         \
                                initial_velocity_z,                                         \
                                initial_inertial_velocity_x,                                \
                                initial_inertial_velocity_y,                                \
                                initial_inertial_velocity_z,                                \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 );",                                      \
                     database )

data1.columns = [ 'rotFrame_vx',                                                            \
                  'rotFrame_vy',                                                            \
                  'rotFrame_vz',                                                            \
                  'inertial_vx',                                                            \
                  'inertial_vy',                                                            \
                  'inertial_vz',                                                            \
                  'launch_azimuth' ]

escape_rotFrame_vx                     = data1[ 'rotFrame_vx' ]
escape_rotFrame_vy                     = data1[ 'rotFrame_vy' ]
escape_rotFrame_vz                     = data1[ 'rotFrame_vz' ]
escape_inertial_vx                     = data1[ 'inertial_vx' ]
escape_inertial_vy                     = data1[ 'inertial_vy' ]
escape_inertial_vz                     = data1[ 'inertial_vz' ]
escape_azimuth                         = data1[ 'launch_azimuth' ]

## get data for re-impact cases
# data2 = pd.read_sql( "SELECT    initial_velocity_x,                                         \
#                                 initial_velocity_y,                                         \
#                                 initial_velocity_z,                                         \
#                                 initial_inertial_velocity_x,                                \
#                                 initial_inertial_velocity_y,                                \
#                                 initial_inertial_velocity_z,                                \
#                                 ROUND( launch_azimuth )                                     \
#                      FROM       regolith_trajectory_results                                 \
#                      WHERE      ( crash_flag = 1 );",                                       \
#                      database )

# data2.columns = [ 'rotFrame_vx',                                                            \
#                   'rotFrame_vy',                                                            \
#                   'rotFrame_vz',                                                            \
#                   'inertial_vx',                                                            \
#                   'inertial_vy',                                                            \
#                   'inertial_vz',                                                            \
#                   'launch_azimuth' ]

# crash_rotFrame_vx                     = data2[ 'rotFrame_vx' ]
# crash_rotFrame_vy                     = data2[ 'rotFrame_vy' ]
# crash_rotFrame_vz                     = data2[ 'rotFrame_vz' ]
# crash_inertial_vx                     = data2[ 'inertial_vx' ]
# crash_inertial_vy                     = data2[ 'inertial_vy' ]
# crash_inertial_vz                     = data2[ 'inertial_vz' ]
# crash_azimuth                         = data2[ 'launch_azimuth' ]

# ## get data for temporary capture cases
# data3 = pd.read_sql( "SELECT    initial_velocity_x,                                         \
#                                 initial_velocity_y,                                         \
#                                 initial_velocity_z,                                         \
#                                 initial_inertial_velocity_x,                                \
#                                 initial_inertial_velocity_y,                                \
#                                 initial_inertial_velocity_z,                                \
#                                 ROUND( launch_azimuth )                                     \
#                      FROM       regolith_trajectory_results                                 \
#                      WHERE      ( end_flag = 1 );",                                         \
#                      database )

# data3.columns = [ 'rotFrame_vx',                                                            \
#                   'rotFrame_vy',                                                            \
#                   'rotFrame_vz',                                                            \
#                   'inertial_vx',                                                            \
#                   'inertial_vy',                                                            \
#                   'inertial_vz',                                                            \
#                   'launch_azimuth' ]

# capture_rotFrame_vx                     = data3[ 'rotFrame_vx' ]
# capture_rotFrame_vy                     = data3[ 'rotFrame_vy' ]
# capture_rotFrame_vz                     = data3[ 'rotFrame_vz' ]
# capture_inertial_vx                     = data3[ 'inertial_vx' ]
# capture_inertial_vy                     = data3[ 'inertial_vy' ]
# capture_inertial_vz                     = data3[ 'inertial_vz' ]
# capture_azimuth                         = data3[ 'launch_azimuth' ]

## plot launch azimuth versus rot frame initial velocity for escape cases
escape_rotFrameInitialVelocity = np.sqrt( escape_rotFrame_vx**2 + escape_rotFrame_vy**2 + escape_rotFrame_vz**2 )
ax1Handle1 = ax1.scatter( escape_azimuth, escape_rotFrameInitialVelocity, color=colors.cnames['purple'], label='actual speed' )

## Plot launch azimuth versus inertial initial velocity for escape cases
escape_inertialInitialVelocity = np.sqrt( escape_inertial_vx**2 + escape_inertial_vy**2 + escape_inertial_vz**2 )
ax2Handle1 = ax2.scatter( escape_azimuth, escape_inertialInitialVelocity, color=colors.cnames['purple'], label='actual speed' )

ax1Handle2, = ax1.plot( directionalEscapeAzimuth, directionalEscapeSpeed, color=colors.cnames['red'], label='directional escape speed', lw=2 )
ax2Handle2, = ax2.plot( directionalEscapeAzimuth, inertialDirectionalEscapeSpeed, color=colors.cnames['red'], label='directional escape speed', lw=2 )

## format axis and title
ax1.set_xlabel('Launch azimuth [deg]')
ax1.set_ylabel('$V_{initial}$ [m/s]')
ax1.set_xlim( 0, 360 )
ax1.set_title( 'Regolith (escape) initial velocity versus launch azimuth (rotating frame)' )
ax1.legend( )
# ax1.legend( loc='center left', bbox_to_anchor=( 0.995,0.5 ) )
ax1.grid( )

## format axis and title
ax2.set_xlabel('Launch azimuth [deg]')
ax2.set_ylabel('$V_{initial}$ [m/s]')
ax2.set_xlim( 0, 360 )
ax2.set_title( 'Regolith (escape) initial velocity versus launch azimuth (inertial frame)' )
ax2.legend( )
# ax2.legend( loc='center left', bbox_to_anchor=( 0.995,0.5 ) )
ax2.grid( )

## Show the plot
# plt.legend( handles = [ ax2Handle1, ax2Handle2, escapeBoundHandle, crashBoundHandle ],
#             bbox_to_anchor = ( 0.95, 0.5 ),
#             loc = 1 )
# plt.legend( bbox_to_anchor=(1.0, 1), loc=1, borderaxespad=0. )
# plt.tight_layout( )
# plt.grid( )
plt.show( )
# plt.savefig('../data/trajectory_for_different_launch_azimuth/trajectoryPlot_'+str(figureIndex+1)+'.png', dpi=300)

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
