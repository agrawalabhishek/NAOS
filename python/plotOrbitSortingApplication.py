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
import mayavi.mlab as mLab
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

def extractDataFromSQL( databaseConnect,
                        currentSolarPhase ):
    # general function for data extraction
    data = pd.read_sql( "SELECT    trajectory_id,                                              \
                                   initial_position_x,                                         \
                                   initial_position_y,                                         \
                                   initial_position_z,                                         \
                                   position_x,                                                 \
                                   position_y,                                                 \
                                   position_z,                                                 \
                                   ROUND( initial_velocity_magnitude ),                        \
                                   ROUND( launch_azimuth )                                     \
                        FROM       regolith_trajectory_results                                 \
                        WHERE      time > 86400.0                                              \
                        AND        ROUND( initial_solar_phase_angle ) = "
                                   + str( currentSolarPhase ) + "                              \
                        ORDER BY   ROUND( initial_velocity_magnitude ) ASC;",                  \
                        databaseConnect )

    data.columns = [ 'trajectory_id',                                                          \
                     'init_pos_x',                                                             \
                     'init_pos_y',                                                             \
                     'init_pos_z',                                                             \
                     'x',                                                                      \
                     'y',                                                                      \
                     'z',                                                                      \
                     'init_vel_mag',                                                           \
                     'launch_azimuth' ]

    x_data                   = data[ 'x' ]
    y_data                   = data[ 'y' ]
    z_data                   = data[ 'z' ]
    initial_velocity_data    = data[ 'init_vel_mag' ]
    azimuth_data             = data[ 'launch_azimuth' ]

    xPositionStart_data      = data[ 'init_pos_x' ]
    yPositionStart_data      = data[ 'init_pos_y' ]
    zPositionStart_data      = data[ 'init_pos_z' ]

    trajectoryID             = data[ 'trajectory_id' ]

    return trajectoryID, x_data, y_data, z_data, initial_velocity_data, azimuth_data,                       \
           xPositionStart_data, yPositionStart_data, zPositionStart_data

def crashMapScatterPlot( unique_reimpact_velocities,
                         initial_velocity,
                         endLatitude,
                         endLongitude,
                         colormap,
                         markertype,
                         plotLabel ):
    for index in range( 0, len( unique_reimpact_velocities ) ):
        current_velocity = unique_reimpact_velocities[ index ]
        current_velocity_indices = np.where( initial_velocity == current_velocity )
        current_velocity_indices = current_velocity_indices[ 0 ]
        plotLatitude = endLatitude[ current_velocity_indices ]
        plotLongitude = endLongitude[ current_velocity_indices ]

        ax1.scatter( plotLongitude, plotLatitude, s=5, c=colormap[index],             \
                     edgecolors='face',                                               \
                     marker=markertype,                                               \
                     label=plotLabel )

def areaToMassRatio( radius, density ):
    crossSectionalArea = np.pi * radius**2
    mass = ( 4.0 / 3.0 ) * np.pi * radius**3 * density
    ratio = crossSectionalArea / mass
    return ratio

# Start timer.
start_time = time.time( )

## Constants for the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0

commonPath = "/media/abhishek/Ashish/Thesis_Simulation_Databases/trailing_edge_perturbations_CDE/"

try:
    database1 = sqlite3.connect(commonPath + "3.2Density_1cmRadius/trailingEdge_3P2Density_1cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database4 = sqlite3.connect(commonPath + "7.5Density_5cmRadius/trailingEdge_7P5Density_5cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

## Operations
solarPhase = [ 45.0, 135.0, 225.0, 315.0 ]
# solarPhase = [ 45.0 ]
fig = plt.figure( figsize=( 15, 15 ) )
gs = gridspec.GridSpec( 2, 2 )
ax2 = plt.subplot( gs[ 0 ] )
ax3 = plt.subplot( gs[ 1 ] )
ax4 = plt.subplot( gs[ 2 ] )
ax5 = plt.subplot( gs[ 3 ] )

plotHandles = [ax2, ax3, ax4, ax5]
customcolormap = plt.cm.Vega20( np.linspace( 0, 1, 16 ) )

for index in range( 0, len(solarPhase) ):
    currentSolarPhase = solarPhase[ index ]
    print "processing data for solar phase = " + str(currentSolarPhase) + " [deg]..."
    ax1 = plotHandles[ index ]


    ( trajID_data1, x_data1, y_data1, z_data1, initial_velocity_data1, azimuth_data1,
        xPositionStart_data1, yPositionStart_data1, zPositionStart_data1 ) = extractDataFromSQL(
                                                                                database1,
                                                                                currentSolarPhase )
    ratio = areaToMassRatio( 1.0e-2, 3.2e3 )
    label1 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    uniqueTrajID = np.unique( trajID_data1 )

    for trajIndex in range( 0, len( uniqueTrajID ) ):
        currentTrajID = uniqueTrajID[ trajIndex ]
        dataIndices = np.where( trajID_data1 == currentTrajID )
        dataIndices = dataIndices[ 0 ]

        currentAzimuth = azimuth_data1[ dataIndices ].tolist( )
        currentVelocity = initial_velocity_data1[ dataIndices ].tolist( )

        ax1.scatter( currentAzimuth[ 0 ], currentVelocity[ 0 ],                                 \
                     s=8, c='blue',                                                             \
                     edgecolors='face',                                                         \
                     marker='^',                                                                \
                     label=label1 if trajIndex == 0 else "__nolegend__" )


    ( trajID_data4, x_data4, y_data4, z_data4, initial_velocity_data4, azimuth_data4,
        xPositionStart_data4, yPositionStart_data4, zPositionStart_data4 ) = extractDataFromSQL(
                                                                                database4,
                                                                                currentSolarPhase )
    ratio = areaToMassRatio( 5.0e-2, 7.5e3 )
    label4 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    uniqueTrajID = np.unique( trajID_data4 )

    for trajIndex in range( 0, len( uniqueTrajID ) ):
        currentTrajID = uniqueTrajID[ trajIndex ]
        dataIndices = np.where( trajID_data4 == currentTrajID )
        dataIndices = dataIndices[ 0 ]

        currentAzimuth = azimuth_data4[ dataIndices ].tolist( )
        currentVelocity = initial_velocity_data4[ dataIndices ].tolist( )

        ax1.scatter( currentAzimuth[ 0 ], currentVelocity[ 0 ],                                      \
                     s=5, c=colors.cnames['red'],                                                    \
                     edgecolors='face',                                                              \
                     marker='o',                                                                     \
                     label=label4 if trajIndex == 0 else "__nolegend__" )

    ax1.set_xlabel('Launch azimuth [deg]')
    ax1.set_ylabel('Launch velocity [m/s]')
    ax1.set_title("Solar phase = " + str(currentSolarPhase) + " [deg]")
    ax1.grid( True )
    ax1.legend( markerscale=3 ).draggable( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
# plt.suptitle( 'Regolith crash map for different regolith types, Ellipsoid trailing edge' )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
