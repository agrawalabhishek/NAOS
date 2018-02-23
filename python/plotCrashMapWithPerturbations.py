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
                               + "multiple_launch_velocity_with_perturbations/"
                               + "simulation_time_9_months/"
                               + "longestEdge_3P2Density_1cmRadius.db")
                               # + "3.2Density_0.4cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

data = pd.read_sql( "SELECT    ROUND( initial_solar_phase_angle )                          \
                    FROM       regolith_trajectory_results                                 \
                    WHERE      ( crash_flag = 1 );",                                       \
                    database )

data.columns = [ 'initial_solar_phase_angle' ]

solarPhases = data[ 'initial_solar_phase_angle' ]
uniqueSolarPhases = np.unique( solarPhases )

fig = plt.figure( figsize=(12, 12) )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )

plotHandles = [ ax1, ax2, ax3, ax4 ]

colors = plt.cm.Vega20( np.linspace( 0, 1, 16 ) )

for solarIndex in range( 0, len( uniqueSolarPhases ) ):
    currentPlotHandle = plotHandles[ solarIndex ]
    currentSolarPhase = uniqueSolarPhases[ solarIndex ]

    data1 = pd.read_sql( "SELECT    initial_position_x,                                         \
                                    initial_position_y,                                         \
                                    initial_position_z,                                         \
                                    position_x,                                                 \
                                    position_y,                                                 \
                                    position_z,                                                 \
                                    ROUND( initial_velocity_magnitude ),                        \
                                    ROUND( launch_azimuth )                                     \
                         FROM       regolith_trajectory_results                                 \
                         WHERE      ( crash_flag = 1 )                                          \
                         AND        ROUND( initial_solar_phase_angle ) = "
                                    + str( currentSolarPhase ) + ";",                           \
                         database )

    data1.columns = [ 'init_pos_x',                                          \
                      'init_pos_y',                                          \
                      'init_pos_z',                                          \
                      'x',                                                   \
                      'y',                                                   \
                      'z',                                                   \
                      'init_vel_mag',                                        \
                      'launch_azimuth' ]

    print "Processing data for solar phase angle = " + str( currentSolarPhase ) + " [deg]\n"

    x                   = data1[ 'x' ]
    y                   = data1[ 'y' ]
    z                   = data1[ 'z' ]
    initial_velocity    = data1[ 'init_vel_mag' ]
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

    unique_reimpact_velocities = np.unique( initial_velocity )

    for index in range( 0, len( unique_reimpact_velocities ) ):
        current_velocity = unique_reimpact_velocities[ index ]
        current_velocity_indices = np.where( initial_velocity == current_velocity )
        current_velocity_indices = current_velocity_indices[ 0 ]
        plotLatitude = endLatitude[ current_velocity_indices ]
        plotLongitude = endLongitude[ current_velocity_indices ]

        plotLatitude  = endLatitude[ current_velocity_indices ]
        plotLongitude = endLongitude[ current_velocity_indices ]

        currentPlotHandle.scatter( plotLongitude, plotLatitude, s=5, c=colors[ index ],           \
                                   edgecolors='face',                                             \
                                   label='$V_{launch}$ = ' + str( current_velocity ) + ' [m/s]' )

    currentPlotHandle.legend( markerscale=7 ).draggable( )
    currentPlotHandle.set_xlabel('longitude [degree]')
    currentPlotHandle.set_ylabel('latitude [degree]')
    # cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )
    currentPlotHandle.set_xlim( -180.0, 180.0 )
    currentPlotHandle.set_yticks( np.arange( -90.0, 90.0, 10.0 ) )
    currentPlotHandle.set_title( 'Solar phase angle = ' + str( currentSolarPhase ) + ' [deg]' )
    currentPlotHandle.grid( True )

if database:
    database.close( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.suptitle( 'Regolith crash map for multiple launch velocities \n Ellipsoid longest edge' )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
