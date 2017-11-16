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
                                   ROUND( launch_azimuth )                                     \
                        FROM       regolith_trajectory_results                                 \
                        WHERE      ( crash_flag = 1 )                                          \
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
                     'launch_azimuth' ]

    x_data                   = data[ 'x' ]
    y_data                   = data[ 'y' ]
    z_data                   = data[ 'z' ]
    initial_velocity_data    = data[ 'init_vel_mag' ]
    azimuth_data             = data[ 'launch_azimuth' ]

    xPositionStart_data      = data[ 'init_pos_x' ]
    yPositionStart_data      = data[ 'init_pos_y' ]
    zPositionStart_data      = data[ 'init_pos_z' ]

    if databaseConnect:
        databaseConnect.close( )

    return x_data, y_data, z_data, initial_velocity_data, azimuth_data,                    \
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

try:
    database1 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "3.2Density_0.1cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database2 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "7.5Density_1cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database3 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "3.2Density_5cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database4 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "7.5Density_5cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database5 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "3.2Density_10cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database6 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "7.5Density_10cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

## Operations
currentSolarPhase = 45.0

customcolormap = plt.cm.Vega20( np.linspace( 0, 1, 16 ) )

fig = plt.figure( figsize=( 10, 7 ) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

plotDatabase1 = True
plotDatabase2 = False
plotDatabase3 = False
plotDatabase4 = False
plotDatabase5 = False
plotDatabase6 = True

if plotDatabase1 == True:
    ( x_data1, y_data1, z_data1, initial_velocity_data1, azimuth_data1,
        xPositionStart_data1, yPositionStart_data1, zPositionStart_data1 ) = extractDataFromSQL(
                                                                                database1,
                                                                                currentSolarPhase )
    ratio = areaToMassRatio( 0.1e-2, 3.2e3 )
    # label1 = 'Olivine, $\rho$ = 3.2 [$g/cm^3$], R = 1.0 [cm], A/M = ' + str(ratio) + ' [$m^2/kg$]'
    label1 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    ## calculate lat long for end points for data1
    endRadialDistance_data1 = np.sqrt( x_data1**2 + y_data1**2 + z_data1**2 )
    endLongitude_data1 = np.arctan2( y_data1, x_data1 ) * 180 / np.pi
    endLatitude_data1 = np.arcsin( z_data1 / endRadialDistance_data1 ) * 180 / np.pi

    ax1.scatter( endLongitude_data1, endLatitude_data1, s=5, c='blue',                           \
                 edgecolors='face',                                                              \
                 marker='^',                                                                     \
                 label=label1 )

if plotDatabase2 == True:
    ( x_data2, y_data2, z_data2, initial_velocity_data2, azimuth_data2,
        xPositionStart_data2, yPositionStart_data2, zPositionStart_data2 ) = extractDataFromSQL(
                                                                                database2,
                                                                                currentSolarPhase )
    ratio = areaToMassRatio( 1.0e-2, 7.5e3 )
    label2 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    ## calculate lat long for end points for data2
    endRadialDistance_data2 = np.sqrt( x_data2**2 + y_data2**2 + z_data2**2 )
    endLongitude_data2 = np.arctan2( y_data2, x_data2 ) * 180 / np.pi
    endLatitude_data2 = np.arcsin( z_data2 / endRadialDistance_data2 ) * 180 / np.pi

    ax1.scatter( endLongitude_data2, endLatitude_data2, s=5, c='red',                            \
                 edgecolors='face',                                                              \
                 marker='o',                                                                     \
                 label=label2 )

if plotDatabase3 == True:
    ( x_data3, y_data3, z_data3, initial_velocity_data3, azimuth_data3,
        xPositionStart_data3, yPositionStart_data3, zPositionStart_data3 ) = extractDataFromSQL(
                                                                                database3,
                                                                                currentSolarPhase )
    ratio = areaToMassRatio( 5.0e-2, 3.2e3 )
    label3 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    ## calculate lat long for end points for data3
    endRadialDistance_data3 = np.sqrt( x_data3**2 + y_data3**2 + z_data3**2 )
    endLongitude_data3 = np.arctan2( y_data3, x_data3 ) * 180 / np.pi
    endLatitude_data3 = np.arcsin( z_data3 / endRadialDistance_data3 ) * 180 / np.pi

    ax1.scatter( endLongitude_data3, endLatitude_data3, s=5, c='green',                          \
                 edgecolors='face',                                                              \
                 marker='*',                                                                     \
                 label=label3 )

if plotDatabase4 == True:
    ( x_data4, y_data4, z_data4, initial_velocity_data4, azimuth_data4,
        xPositionStart_data4, yPositionStart_data4, zPositionStart_data4 ) = extractDataFromSQL(
                                                                                database4,
                                                                                currentSolarPhase )
    ratio = areaToMassRatio( 5.0e-2, 7.5e3 )
    label4 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    ## calculate lat long for end points for data4
    endRadialDistance_data4 = np.sqrt( x_data4**2 + y_data4**2 + z_data4**2 )
    endLongitude_data4 = np.arctan2( y_data4, x_data4 ) * 180 / np.pi
    endLatitude_data4 = np.arcsin( z_data4 / endRadialDistance_data4 ) * 180 / np.pi

    ax1.scatter( endLongitude_data4, endLatitude_data4, s=5, c=colors.cnames['darkorange'],      \
                 edgecolors='face',                                                              \
                 marker='+',                                                                     \
                 label=label4 )

if plotDatabase5 == True:
    ( x_data5, y_data5, z_data5, initial_velocity_data5, azimuth_data5,
        xPositionStart_data5, yPositionStart_data5, zPositionStart_data5 ) = extractDataFromSQL(
                                                                                database5,
                                                                                currentSolarPhase )
    ratio = areaToMassRatio( 10e-2, 3.2e3 )
    label5 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    # calculate lat long for end points for data5
    endRadialDistance_data5 = np.sqrt( x_data5**2 + y_data5**2 + z_data5**2 )
    endLongitude_data5 = np.arctan2( y_data5, x_data5 ) * 180 / np.pi
    endLatitude_data5 = np.arcsin( z_data5 / endRadialDistance_data5 ) * 180 / np.pi

    ax1.scatter( endLongitude_data5, endLatitude_data5, s=5, c=colors.cnames['magenta'],         \
                 edgecolors='face',                                                              \
                 marker='s',                                                                     \
                 label=label5 )

if plotDatabase6 == True:
    ( x_data6, y_data6, z_data6, initial_velocity_data6, azimuth_data6,
        xPositionStart_data6, yPositionStart_data6, zPositionStart_data6 ) = extractDataFromSQL(
                                                                                database6,
                                                                                currentSolarPhase )
    ratio = areaToMassRatio( 10e-2, 7.5e3 )
    label6 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

    ## calculate lat long for end points for data6
    endRadialDistance_data6 = np.sqrt( x_data6**2 + y_data6**2 + z_data6**2 )
    endLongitude_data6 = np.arctan2( y_data6, x_data6 ) * 180 / np.pi
    endLatitude_data6 = np.arcsin( z_data6 / endRadialDistance_data6 ) * 180 / np.pi

    ax1.scatter( endLongitude_data6, endLatitude_data6, s=5, c=colors.cnames['purple'],          \
                 edgecolors='face',                                                              \
                 marker='D',                                                                     \
                 label=label6 )

## calculate the lat long for starting point (same for all regoliths type in this case)
# startRadialDistance = np.sqrt( xPositionStart_data1**2 +
#                                yPositionStart_data1**2 +
#                                zPositionStart_data1**2 )
# startLongitude = np.arctan2( yPositionStart_data1, xPositionStart_data1 ) * 180.0 / np.pi
# startLatitude = np.arcsin( zPositionStart_data1 / startRadialDistance ) * 180.0 / np.pi

# unique_reimpact_velocities_data1 = np.unique( initial_velocity_data1 )
# unique_reimpact_velocities_data2 = np.unique( initial_velocity_data2 )
# unique_reimpact_velocities_data3 = np.unique( initial_velocity_data3 )
# unique_reimpact_velocities_data4 = np.unique( initial_velocity_data4 )
# unique_reimpact_velocities_data5 = np.unique( initial_velocity_data5 )
# unique_reimpact_velocities_data6 = np.unique( initial_velocity_data6 )

# crashMapScatterPlot( unique_reimpact_velocities_data1,
#                      initial_velocity_data1,
#                      endLatitude_data1,
#                      endLongitude_data1,
#                      colors,
#                      markertype='o',
#                      plotLabel=label1 )

# crashMapScatterPlot( unique_reimpact_velocities_data2,
#                      initial_velocity_data2,
#                      endLatitude_data2,
#                      endLongitude_data2,
#                      colors,
#                      markertype='^',
#                      plotLabel=label2 )

# crashMapScatterPlot( unique_reimpact_velocities_data3,
#                      initial_velocity_data3,
#                      endLatitude_data3,
#                      endLongitude_data3,
#                      colors,
#                      markertype='*',
#                      plotLabel=label3 )

############################ Experiment - contours of initial launch velocity ######################

# sort data based on velocity
# sortedLongitude_data1 = []
# sortedLatitude_data1 = []
# sortedVelocity_data1 = []
# for index in range( 0, len( unique_reimpact_velocities_data1 ) ):
#         current_velocity = unique_reimpact_velocities_data1[ index ]
#         current_velocity_indices = np.where( initial_velocity_data1 == current_velocity )
#         current_velocity_indices = current_velocity_indices[ 0 ]
#         sortedLatitude_data1.append( endLatitude_data1[ current_velocity_indices ].tolist() )
#         sortedLongitude_data1.append( endLongitude_data1[ current_velocity_indices ].tolist() )
#         sortedVelocity_data1.append( initial_velocity_data1[ current_velocity_indices ].tolist() )

# xData = sum( sortedLongitude_data1, [ ] )
# yData = sum( sortedLatitude_data1, [ ] )
# zData = sum( sortedVelocity_data1, [ ] )

# xData = endLongitude_data1
# yData = endLatitude_data1
# zData = initial_velocity_data1

# xi = np.linspace( -180.0, 180.0, 1000 )
# yi = np.linspace( -90.0, 90.0, 1000 )
# zi = mlab.griddata( xData, yData, zData, xi, yi, interp='linear' )
# ax1.contour( xi, yi, zi, 16,
#              linewidths=0.5, colors='black', label='Launch velocity contour lines' )

# grid_x, grid_y = np.mgrid[ -180:180:100j, -90:90:100j ]
# zi = griddata( (xData, yData), zData, (grid_x, grid_y), method='cubic' )
# contourplot = ax1.contourf( grid_x.T, grid_y.T, zi.T,
#                            16,
#                            linewidths=1.0,
#                            cmap=plt.cm.get_cmap("winter"),
#                            label='Launch velocity contour lines' )
# fig.colorbar( contourplot, ax=ax1 )

ax1.set_xlabel('longitude [degree]')
ax1.set_ylabel('latitude [degree]')
# cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )
ax1.set_xlim( -180.0, 180.0 )
ax1.set_yticks( np.arange( -90.0, 90.0, 15.0 ) )
# ax1.set_title( 'Solar phase angle = ' + str( currentSolarPhase ) + ' [deg]' )
ax1.grid( True )
ax1.legend( markerscale=2 ).draggable( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.suptitle( 'Regolith crash map for different regolith types \n Ellipsoid longest edge' +
              '\n Solar phase = ' + str(currentSolarPhase) + ' [deg]' )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
