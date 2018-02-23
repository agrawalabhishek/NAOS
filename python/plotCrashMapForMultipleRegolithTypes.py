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

## Constants for the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0

## Operations
uniqueSolarPhases = range( 45, 360, 90 )

fig = plt.figure( figsize=( 12, 12 ) )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )

plotHandles = [ ax1, ax2, ax3, ax4 ]

colors = plt.cm.Vega20( np.linspace( 0, 1, 16 ) )

try:
    database1 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "3.2Density_0.4cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database2 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "7.5Density_4cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

for solarIndex in range( 0, len( uniqueSolarPhases ) ):
    currentPlotHandle = plotHandles[ solarIndex ]
    currentSolarPhase = uniqueSolarPhases[ solarIndex ]

    print "Extracting data now...\n"

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
                         database1 )

    data1.columns = [ 'init_pos_x',                                          \
                      'init_pos_y',                                          \
                      'init_pos_z',                                          \
                      'x',                                                   \
                      'y',                                                   \
                      'z',                                                   \
                      'init_vel_mag',                                        \
                      'launch_azimuth' ]

    x_data1                   = data1[ 'x' ]
    y_data1                   = data1[ 'y' ]
    z_data1                   = data1[ 'z' ]
    initial_velocity_data1    = data1[ 'init_vel_mag' ]
    azimuth_data1             = data1[ 'launch_azimuth' ]

    xPositionStart_data1      = data1[ 'init_pos_x' ]
    yPositionStart_data1      = data1[ 'init_pos_y' ]
    zPositionStart_data1      = data1[ 'init_pos_z' ]

    data2 = pd.read_sql( "SELECT    initial_position_x,                                         \
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
                         database2 )

    data2.columns = [ 'init_pos_x',                                          \
                      'init_pos_y',                                          \
                      'init_pos_z',                                          \
                      'x',                                                   \
                      'y',                                                   \
                      'z',                                                   \
                      'init_vel_mag',                                        \
                      'launch_azimuth' ]

    x_data2                   = data2[ 'x' ]
    y_data2                   = data2[ 'y' ]
    z_data2                   = data2[ 'z' ]
    initial_velocity_data2    = data2[ 'init_vel_mag' ]
    azimuth_data2             = data2[ 'launch_azimuth' ]

    xPositionStart_data2      = data2[ 'init_pos_x' ]
    yPositionStart_data2      = data2[ 'init_pos_y' ]
    zPositionStart_data2      = data2[ 'init_pos_z' ]

    print "Processing data for solar phase angle = " + str( currentSolarPhase ) + " [deg]\n"

    ## calculate the lat long for starting point (same for all regoliths type in this case)
    startRadialDistance = np.sqrt( xPositionStart_data1**2 +
                                   yPositionStart_data1**2 +
                                   zPositionStart_data1**2 )
    startLongitude = np.arctan2( yPositionStart_data1, xPositionStart_data1 ) * 180.0 / np.pi
    startLatitude = np.arcsin( zPositionStart_data1 / startRadialDistance ) * 180.0 / np.pi

    ## calculate lat long for end points for data1
    endRadialDistance_data1 = np.sqrt( x_data1**2 + y_data1**2 + z_data1**2 )
    endLongitude_data1 = np.arctan2( y_data1, x_data1 ) * 180 / np.pi
    endLatitude_data1 = np.arcsin( z_data1 / endRadialDistance_data1 ) * 180 / np.pi

    ## calculate lat long for end points for data2
    endRadialDistance_data2 = np.sqrt( x_data2**2 + y_data2**2 + z_data2**2 )
    endLongitude_data2 = np.arctan2( y_data2, x_data2 ) * 180 / np.pi
    endLatitude_data2 = np.arcsin( z_data2 / endRadialDistance_data2 ) * 180 / np.pi

    unique_reimpact_velocities_data1 = np.unique( initial_velocity_data1 )
    unique_reimpact_velocities_data2 = np.unique( initial_velocity_data2 )


    currentPlotHandle.scatter( endLongitude_data1, endLatitude_data1, s=5, c='black',           \
                               edgecolors='face',                                               \
                               marker='o',                                                      \
                               label='Olivine, Density = 3.2 [$g/cm^3$], Radius = 1.0 [cm]' )

    currentPlotHandle.scatter( endLongitude_data2, endLatitude_data2, s=5, c='red',             \
                               edgecolors='face',                                               \
                               marker='^',                                                      \
                               label='Fe-Ni, Density = 7.5 [$g/cm^3$], Radius = 1.0 [cm]' )

    currentPlotHandle.set_xlabel('longitude [degree]')
    currentPlotHandle.set_ylabel('latitude [degree]')
    # cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )
    currentPlotHandle.set_xlim( -180.0, 180.0 )
    currentPlotHandle.set_yticks( np.arange( -90.0, 90.0, 10.0 ) )
    currentPlotHandle.set_title( 'Solar phase angle = ' + str( currentSolarPhase ) + ' [deg]' )
    currentPlotHandle.grid( True )

currentPlotHandle.legend( markerscale=2 ).draggable( )

if database1:
    database1.close( )
if database2:
    database2.close( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.suptitle( 'Regolith crash map for different regolith types \n Ellipsoid longest edge' )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
