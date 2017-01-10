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
        database = sqlite3.connect("../data/regolith_launched_from_leading_edge/single_launch_velocity/leadingEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

data1 = pd.read_sql( "SELECT    position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                velocity_x,                                                 \
                                velocity_y,                                                 \
                                velocity_z,                                                 \
                                time,                                                       \
                                escape_flag,                                                \
                                sma,                                                        \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 );",                                      \
                     database )

data1.columns = [ 'x',                                                   \
                  'y',                                                   \
                  'z',                                                   \
                  'vx',                                                  \
                  'vy',                                                  \
                  'vz',                                                  \
                  'time',                                                \
                  'escape_flag',                                         \
                  'sma',                                                 \
                  'launch_azimuth' ]

escape_x           = data1[ 'x' ]
escape_y           = data1[ 'y' ]
escape_z           = data1[ 'z' ]
escape_vx          = data1[ 'vx' ]
escape_vy          = data1[ 'vy' ]
escape_vz          = data1[ 'vz' ]
escape_t           = data1[ 'time' ]
escape_sma         = data1[ 'sma' ]
escape_azimuth     = data1[ 'launch_azimuth' ]


# get data for reimpact cases
data2 = pd.read_sql( "SELECT    position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                velocity_x,                                                 \
                                velocity_y,                                                 \
                                velocity_z,                                                 \
                                time,                                                       \
                                crash_flag,                                                 \
                                sma,                                                        \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( crash_flag = 1 );",                                      \
                     database )

data2.columns = [ 'x',                                                   \
                  'y',                                                   \
                  'z',                                                   \
                  'vx',                                                  \
                  'vy',                                                  \
                  'vz',                                                  \
                  'time',                                                \
                  'crash_flag',                                          \
                  'sma',                                                 \
                  'launch_azimuth' ]

crash_x           = data2[ 'x' ]
crash_y           = data2[ 'y' ]
crash_z           = data2[ 'z' ]
crash_vx          = data2[ 'vx' ]
crash_vy          = data2[ 'vy' ]
crash_vz          = data2[ 'vz' ]
crash_t           = data2[ 'time' ]
crash_sma         = data2[ 'sma' ]
crash_azimuth     = data2[ 'launch_azimuth' ]

## plot launch azimuth versus crash velocity
crashVelocity = np.sqrt( crash_vx**2 + crash_vy**2 + crash_vz**2 )
ax1.scatter( crash_azimuth, crashVelocity, color=colors.cnames['green'] )

## format axis and title
ax1.set_xlabel('Launch azimuth [deg]')
ax1.set_ylabel('$V_{impact}$ [m/s]')
ax1.set_xlim( 0, 360 )
ax1.set_title( 'Regolith launch azimuth versus $V_{impact}$' )
# ax1.legend( loc='center left', bbox_to_anchor=( 0.995,0.5 ) )
ax1.grid( )

## Plot launch azimuth versus hev
# get hev first
hev = np.sqrt( mu / ( -1.0 * escape_sma ) )
ax2.scatter( escape_azimuth, hev, color=colors.cnames['purple'] )

## format axis and title
ax2.set_xlabel('Launch azimuth [deg]')
ax2.set_ylabel('$V_\infty$ [m/s]')
ax2.set_xlim( 0, 360 )
ax2.set_title( 'Regolith launch azimuth versus $V_\infty$' )
# ax2.legend( loc='center left', bbox_to_anchor=( 0.995,0.5 ) )
ax2.grid( )

## Show the plot
# plt.legend( )
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
