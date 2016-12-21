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

## Set up the figure
fig = plt.figure( )
gs = gridspec.GridSpec( 2, 1, height_ratios = [ 2.5, 1 ] )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
# ax1 = fig.add_subplot(111)
# ax2 = fig.add_subplot(212)

## Operations
# Connect to SQLite database.
try:
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

# for angleValue in tqdm(angleArray):
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

x           = data1[ 'x' ]
y           = data1[ 'y' ]
z           = data1[ 'z' ]
vx          = data1[ 'vx' ]
vy          = data1[ 'vy' ]
vz          = data1[ 'vz' ]
t           = data1[ 'time' ]
escapeFlag  = data1[ 'escape_flag' ]
sma         = data1[ 'sma' ]
azimuth     = data1[ 'launch_azimuth' ]

# all ejecta have the same starting point in this particular case
data2 = pd.read_sql( "SELECT     position_x,                                                 \
                                 position_y,                                                 \
                                 position_z,                                                 \
                                 velocity_x,                                                 \
                                 velocity_y,                                                 \
                                 velocity_z,                                                 \
                                 time,                                                       \
                                 ROUND( launch_azimuth )                                     \
                     FROM        regolith_trajectory_results                                 \
                     WHERE       ( start_flag = 1 )                                          \
                     AND         ( ROUND( launch_azimuth ) = 0 );",                          \
                     database )

data2.columns = [ 'x',                                                   \
                  'y',                                                   \
                  'z',                                                   \
                  'vx',                                                  \
                  'vy',                                                  \
                  'vz',                                                  \
                  'time',                                                \
                  'launch_azimuth' ]

xPositionStart = data2[ 'x' ]
yPositionStart = data2[ 'y' ]
zPositionStart = data2[ 'z' ]
xVelocityStart = data2[ 'vx' ]
yVelocityStart = data2[ 'vy' ]
zVelocityStart = data2[ 'vz' ]

# get the unique declination angle
data3 = pd.read_sql( "SELECT     DISTINCT launch_declination                                 \
                     FROM        regolith_trajectory_results                                 \
                     WHERE       ( crash_flag = 1 );",                                       \
                     database )

data3.columns = [ 'declination' ]
launchDeclination = data3[ 'declination' ]

asteroidPeriod = 2 * np.pi / Wz

## Plot histogram - number of particles against $V_\infty$
# get the hyperbolic excess velocity (hev)
mu = 876514
hev = np.sqrt( mu / ( -1.0 * sma ) )
ax1.hist( hev, 50, facecolor='green' )

# plotting only until certain 'hev' values
# newHEV = []
# for temp in hev:
#     if temp <= 20: newHEV.append( temp )
# ax1.hist( newHEV, 50, facecolor='green' )

## format axis and title
ax1.set_ylabel('Number of particles')
ax1.set_xlabel('$V_\infty$ [m/s]')
ax1.set_title( 'Histogram - Number of particles v/s $V_\infty$' )
# ax1.legend( loc='center left', bbox_to_anchor=( 0.995,0.5 ) )
ax1.grid( )

## Plot the metadata (launch initial conditions)
ax2.axis( 'off' )
metadata_table = []
metadata_table.append( [ "Initial X coordinate", xPositionStart[0], "[m]" ] )
metadata_table.append( [ "Initial Y coordinate", yPositionStart[0], "[m]" ] )
metadata_table.append( [ "Initial Z coordinate", zPositionStart[0], "[m]" ] )
metadata_table.append( [ "Initial X velocity", xVelocityStart[0], "[m/s]" ] )
metadata_table.append( [ "Initial Y velocity", yVelocityStart[0], "[m/s]" ] )
metadata_table.append( [ "Initial Z velocity", zVelocityStart[0], "[m/s]" ] )
metadata_table.append( [ "Launch declination", launchDeclination[0], "[deg]" ] )
metadata_table.append( [ "Asteroid rotation period", asteroidPeriod, "[s]" ] )
table = ax2.table( cellText = metadata_table, colLabels = None, cellLoc = 'center', loc = 'center' )
table_properties = table.properties( )
table_cells = table_properties[ 'child_artists' ]
for cell in table_cells: cell.set_height( 0.15 )
cell_dict = table.get_celld( )
for row in xrange( 0, 8 ): cell_dict[ ( row, 2 ) ].set_width( 0.1 )

## Show the plot
# plt.legend( )
# plt.tight_layout( )
plt.grid( )
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
