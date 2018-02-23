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
    data = pd.read_sql( "SELECT    initial_position_x,                                         \
                                   initial_position_y,                                         \
                                   initial_position_z,                                         \
                                   position_x,                                                 \
                                   position_y,                                                 \
                                   position_z,                                                 \
                                   ROUND( initial_velocity_magnitude ),                        \
                                   sma,                                                        \
                                   ROUND( launch_azimuth )                                     \
                        FROM       regolith_trajectory_results                                 \
                        WHERE      ( escape_flag = 1 )                                         \
                        AND        ROUND( initial_solar_phase_angle ) = "
                                   + str( currentSolarPhase ) + "                              \
                        ORDER BY   ROUND( initial_velocity_magnitude ) ASC;",                  \
                        databaseConnect )

    data.columns = [ 'init_pos_x',                                                             \
                     'init_pos_y',                                                             \
                     'init_pos_z',                                                             \
                     'x',                                                                      \
                     'y',                                                                      \
                     'z',                                                                      \
                     'init_vel_mag',                                                           \
                     'sma',                                                                    \
                     'launch_azimuth' ]

    x_data                   = data[ 'x' ]
    y_data                   = data[ 'y' ]
    z_data                   = data[ 'z' ]
    initial_velocity_data    = data[ 'init_vel_mag' ]
    azimuth_data             = data[ 'launch_azimuth' ]

    xPositionStart_data      = data[ 'init_pos_x' ]
    yPositionStart_data      = data[ 'init_pos_y' ]
    zPositionStart_data      = data[ 'init_pos_z' ]

    sma = data['sma']

    # if databaseConnect:
    #     databaseConnect.close( )

    return sma, x_data, y_data, z_data, initial_velocity_data, azimuth_data,                    \
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

commonPath = "/media/abhishek/Ashish/Thesis_Simulation_Databases/trailing_edge_perturbations_CDE/"
commonPath2 = "../data/regolith_launched_from_trailing_edge/multiple_launch_velocity_with_perturbations/simulation_time_9_months/"
## Constants for the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0

try:
    database1 = sqlite3.connect(commonPath + "3.2Density_1cmRadius/trailingEdge_3P2Density_1cmRadius.db")
    # database1 = sqlite3.connect(commonPath2 + "/trailingEdge_3P2Density_1cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database2 = sqlite3.connect(commonPath + "7.5Density_1cmRadius/trailingEdge_7P5Density_1cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database3 = sqlite3.connect(commonPath + "3.2Density_5cmRadius/trailingEdge_3P2Density_5cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database4 = sqlite3.connect(commonPath + "7.5Density_5cmRadius/trailingEdge_7P5Density_5cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

## Operations
muAsteroid = 876514
# solarPhase = [ 45.0, 135.0, 225.0, 315.0 ]
solarPhase = [225.0]

customcolormap = plt.cm.Vega20( np.linspace( 0, 1, 16 ) )

for index in range( 0, len(solarPhase) ):
    currentSolarPhase = solarPhase[ index ]
    print "processing data for solar phase = " + str(currentSolarPhase) + " [deg]..."

    fig = plt.figure( figsize=( 15, 15 ) )
    gs = gridspec.GridSpec( 2, 2 )
    ax2 = plt.subplot( gs[ 0 ] )
    ax3 = plt.subplot( gs[ 1 ] )
    ax4 = plt.subplot( gs[ 2 ] )

    ratio = areaToMassRatio( 1.0e-2, 3.2e3 )
    label1 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    ( sma1, x_data1, y_data1, z_data1, initial_velocity_data1, azimuth_data1,
        xPositionStart_data1, yPositionStart_data1, zPositionStart_data1 ) = extractDataFromSQL(
                                                                                database1,
                                                                                currentSolarPhase )

    data_indices = np.where( initial_velocity_data1 == 9.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma1[ data_indices ] )
    ax2.scatter( azimuth_data1[data_indices], hev )

    data_indices = np.where( initial_velocity_data1 == 13.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma1[ data_indices ] )
    ax3.scatter( azimuth_data1[data_indices], hev )

    data_indices = np.where( initial_velocity_data1 == 15.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma1[ data_indices ] )
    ax4.scatter( azimuth_data1[data_indices], hev, label=label1 )

    ratio = areaToMassRatio( 1.0e-2, 7.5e3 )
    label2 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    ( sma2, x_data2, y_data2, z_data2, initial_velocity_data2, azimuth_data2,
        xPositionStart_data2, yPositionStart_data2, zPositionStart_data2 ) = extractDataFromSQL(
                                                                                database2,
                                                                                currentSolarPhase )

    data_indices = np.where( initial_velocity_data2 == 9.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma2[ data_indices ] )
    ax2.scatter( azimuth_data2[data_indices], hev )

    data_indices = np.where( initial_velocity_data2 == 13.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma2[ data_indices ] )
    ax3.scatter( azimuth_data2[data_indices], hev )

    data_indices = np.where( initial_velocity_data2 == 15.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma2[ data_indices ] )
    ax4.scatter( azimuth_data2[data_indices], hev, label=label2 )

    ratio = areaToMassRatio( 5.0e-2, 3.2e3 )
    label3 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    ( sma3, x_data3, y_data3, z_data3, initial_velocity_data3, azimuth_data3,
        xPositionStart_data3, yPositionStart_data3, zPositionStart_data3 ) = extractDataFromSQL(
                                                                                database3,
                                                                                currentSolarPhase )

    data_indices = np.where( initial_velocity_data3 == 9.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma3[ data_indices ] )
    ax2.scatter( azimuth_data3[data_indices], hev )

    data_indices = np.where( initial_velocity_data3 == 13.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma3[ data_indices ] )
    ax3.scatter( azimuth_data3[data_indices], hev )

    data_indices = np.where( initial_velocity_data3 == 15.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma3[ data_indices ] )
    ax4.scatter( azimuth_data3[data_indices], hev, label=label3 )

    ratio = areaToMassRatio( 5.0e-2, 7.5e3 )
    label4 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    ( sma4, x_data4, y_data4, z_data4, initial_velocity_data4, azimuth_data4,
        xPositionStart_data4, yPositionStart_data4, zPositionStart_data4 ) = extractDataFromSQL(
                                                                                database4,
                                                                                currentSolarPhase )

    data_indices = np.where( initial_velocity_data4 == 9.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma4[ data_indices ] )
    ax2.scatter( azimuth_data4[data_indices], hev )

    data_indices = np.where( initial_velocity_data4 == 13.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma4[ data_indices ] )
    ax3.scatter( azimuth_data4[data_indices], hev )

    data_indices = np.where( initial_velocity_data4 == 15.0 )
    data_indices = data_indices[ 0 ]
    hev = np.sqrt( -muAsteroid / sma4[ data_indices ] )
    ax4.scatter( azimuth_data4[data_indices], hev, label=label4 )

    ax2.set_xlabel('Launch azimuth [degree]')
    ax2.set_ylabel('HEV [m/s]')
    ax2.set_title("Launch velocity = 9 [m/s]")
    ax2.grid( True )

    ax3.set_xlabel('Launch azimuth [degree]')
    ax3.set_ylabel('HEV [m/s]')
    ax3.set_title("Launch velocity = 13 [m/s]")
    ax3.grid( True )

    ax4.set_xlabel('Launch azimuth [degree]')
    ax4.set_ylabel('HEV [m/s]')
    ax4.set_title("Launch velocity = 15 [m/s]")
    ax4.grid( True )
    ax4.legend( markerscale=4 ).draggable( )

    plt.suptitle( "Regolith escape behavior with varying launch azimuth\n"
                    + "Solar phase = " + str( currentSolarPhase ) + " [deg]" )


# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
