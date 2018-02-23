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
muAsteroid = 876514
muSun = 1.32712440018 * 1.0e+20
AU = 149597870700.0

def plotCircle( radius, plotHandle ):
        angleRange = np.linspace( 0.0, 2.0*np.pi, 360.0 )
        x = radius * np.cos( angleRange )
        y = radius * np.sin( angleRange )
        plotHandle.plot( x, y )

## Operations
# Connect to SQLite database.
try:
    database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db" )
    # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )
    # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid_with_perturbations/longestEdge.db" )

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

fig = plt.figure( )
gs = gridspec.GridSpec( 1, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )

## get data
#define the angle range for which the data has to be exctracted
launch_azimuth_range = tuple( ( np.arange( 0.0, 360.0, 10.0 ).tolist( ) ) )

data1 = pd.read_sql( "SELECT    trajectory_id,                                                  \
                                crash_flag,                                                     \
                                escape_flag                                                     \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( initial_velocity_magnitude ) = 8                         \
                     AND        ROUND( launch_azimuth ) IN " + str( launch_azimuth_range ) + "  \
                     AND        ( crash_flag = 1 OR escape_flag = 1 );",                        \
                     database )

data1.columns = [ 'traj_id',                                                                    \
                  'crash_flag',                                                                 \
                  'escape_flag' ]

trajectory_id = data1[ 'traj_id' ]
crash_flag    = data1[ 'crash_flag' ]
escape_flag   = data1[ 'escape_flag' ]

trajectory_id_list = trajectory_id.tolist( )
trajectory_id_tuple = tuple( trajectory_id_list )

data2 = pd.read_sql( "SELECT    trajectory_id,                                                  \
                                inclination,                                                    \
                                raan,                                                           \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth )                                         \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      trajectory_id IN " + str( trajectory_id_tuple ) + ";",          \
                     database )

data2.columns = [ 'trajectory_id',                                                      \
                  'inclination',                                                        \
                  'raan',                                                               \
                  'vel_mag',                                                            \
                  'launch_azimuth' ]

data2_trajectory_id = data2[ 'trajectory_id' ]
inclination         = data2[ 'inclination' ]
raan                = data2[ 'raan' ]
velocity_magnitude  = data2[ 'vel_mag' ]
launch_azimuth      = data2[ 'launch_azimuth' ]

if database:
    database.close( )

print "processing data now...\n"

## plot data
unique_launch_azimuth = ( np.unique( launch_azimuth ) )
colors = plt.cm.Vega20( np.linspace( 0, 1, len( unique_launch_azimuth ) ) )
loopCount = 0

for index in ( range( 0, len( unique_launch_azimuth ) ) ):
    loopCount = loopCount + 1
    print "loop Count = " + str( loopCount )
    current_launch_azimuth = unique_launch_azimuth[ index ]
    current_launch_azimuth_indices = np.where( launch_azimuth == current_launch_azimuth )
    current_launch_azimuth_indices = current_launch_azimuth_indices[ 0 ]

    currentTraj_inclination = inclination[ current_launch_azimuth_indices ]
    currentTraj_raan = raan[ current_launch_azimuth_indices ]

    current_trajectory_id = data2_trajectory_id[ current_launch_azimuth_indices ]
    current_trajectory_id = np.unique( current_trajectory_id )
    print "Trajectory ID = " + str( current_trajectory_id[ 0 ] )
    trajectory_id_index = np.where( trajectory_id == current_trajectory_id[ 0 ] )
    trajectory_id_index = trajectory_id_index[ 0 ]

    current_crash_flag_value = crash_flag[ trajectory_id_index ]
    current_crash_flag_value =  current_crash_flag_value.tolist( )
    print "Reimpact flag value = " + str( current_crash_flag_value )

    current_escape_flag_value = escape_flag[ trajectory_id_index ]
    current_escape_flag_value = current_escape_flag_value.tolist( )
    print "Escape flag value = " + str( current_escape_flag_value )
    print "\n"

    i_cos_raan = currentTraj_inclination * np.cos( currentTraj_raan * np.pi / 180.0 )
    i_cos_raan = i_cos_raan.tolist( )
    i_sin_raan = currentTraj_inclination * np.sin( currentTraj_raan * np.pi / 180.0 )
    i_sin_raan = i_sin_raan.tolist( )

    if current_crash_flag_value[ 0 ] == 1:
        # ax1.scatter( i_cos_raan, i_sin_raan, s=5, c=colors[index],    \
        #              edgecolors='face',                         \
        #              label='Launch azimuth = ' + str(current_launch_azimuth) + '[deg]' )
        ax1.scatter( i_cos_raan, i_sin_raan, s=5, c=colors[index],    \
                     edgecolors='face' )
        ax1.scatter( i_cos_raan[ len(i_cos_raan)-1 ], i_sin_raan[ len(i_sin_raan)-1 ],  \
                     s=7, c='black', edgecolors='face' )
        ax1.text( i_cos_raan[0], i_sin_raan[0], str(current_launch_azimuth), \
                  size=10, zorder=1, color='black', rotation='vertical', va='bottom' )
        # ax1.scatter( currentTraj_inclination, currentTraj_raan, s=5, c=colors[index],    \
        #              edgecolors='face',                                                  \
        #              label='Launch azimuth = ' + str(current_launch_azimuth) + '[deg]' )

    elif current_escape_flag_value[ 0 ] == 1:
        # ax2.scatter( i_cos_raan, i_sin_raan, s=5, c=colors[index],    \
        #              edgecolors='face',                         \
        #              label='Launch azimuth = ' + str(current_launch_azimuth) + '[deg]' )
        ax2.scatter( i_cos_raan, i_sin_raan, s=5, c=colors[index],    \
                     edgecolors='face' )
        ax2.scatter( i_cos_raan[ len(i_cos_raan)-1 ], i_sin_raan[ len(i_sin_raan)-1 ], \
                     s=7, c='black', edgecolors='face' )
        ax2.text( i_cos_raan[0], i_sin_raan[0], str(current_launch_azimuth), \
                  size=10, zorder=1, color='black', rotation='vertical', va='bottom' )
        # ax2.scatter( currentTraj_inclination, currentTraj_raan, s=5, c=colors[index],    \
        #              edgecolors='face',                                                  \
        #              label='Launch azimuth = ' + str(current_launch_azimuth) + '[deg]' )

ax1.set_ylabel('$i.sin(\Omega)$')
ax1.set_xlabel('$i.cos(\Omega)$')
ax1.grid(True)
ax1.set_title('Reimpact case')
# ax1.legend( markerscale=10 ).draggable( )

ax2.set_ylabel('$i.sin(\Omega)$')
ax2.set_xlabel('$i.cos(\Omega)$')
ax2.grid(True)
ax2.set_title('Escape case')
# ax2.legend( markerscale=10 ).draggable( )

## close the database
if database:
    database.close( )

## set global plot title
plt.suptitle( 'Inclination vector plot, $V_{launch}$ = ' + str( velocity_magnitude[ 0 ] ) + ' [m/s]' )

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
