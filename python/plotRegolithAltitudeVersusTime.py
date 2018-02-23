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

## Plot the ellipsoidal shape of the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593

# Connect to SQLite database.
try:
        # database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )
        # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db")
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   # + "3.2Density_1cmSize/"
                                   + "longestEdge_3P2Density_1cmRadius.db")
                                   # + "longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

# azimuthRange = tuple( ( np.arange( 0.0, 360.0, 1.0 ) ).tolist( ) )

# data1 = pd.read_sql( "SELECT     trajectory_id                                                   \
#                       FROM       regolith_trajectory_results                                     \
#                       WHERE      ROUND( launch_azimuth ) IN " + str(azimuthRange) + "            \
#                       AND        ROUND( initial_velocity_magnitude ) = 7.0                       \
#                       AND        end_flag = 1;",                                                 \
#                       database )

data1 = pd.read_sql( "SELECT     trajectory_id                                                   \
                      FROM       regolith_trajectory_results                                     \
                      WHERE      end_flag = 1;",                                                 \
                      database )

data1.columns = [ 'trajectory_id' ]
trajectory_id = data1[ 'trajectory_id' ]
trajectory_id = tuple( trajectory_id.tolist( ) )

data = pd.read_sql( "SELECT     trajectory_id,                                                  \
                                position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth ),                                        \
                                ROUND( initial_solar_phase_angle ),                             \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      trajectory_id IN " + str( trajectory_id ) + ";",                \
                     database )

if database:
    database.close( )

data.columns = [ 'trajectory_id',                                                               \
                 'x',                                                                           \
                 'y',                                                                           \
                 'z',                                                                           \
                 'velocity_magnitude',                                                          \
                 'launch_azimuth',                                                              \
                 'initial_solar_phase',                                                         \
                 'time' ]

trajectory_id       = data[ 'trajectory_id' ]
x                   = data[ 'x' ]
y                   = data[ 'y' ]
z                   = data[ 'z' ]
velocityMagnitude   = data[ 'velocity_magnitude' ]
launchAzimuth       = data[ 'launch_azimuth' ]
solarPhase          = data[ 'initial_solar_phase' ]
t                   = data[ 'time' ]

print "Processing data now...\n"

# fig = plt.figure( )
# gs = gridspec.GridSpec( 2, 2 )
# ax1 = plt.subplot( gs[ 0 ] )
# ax2 = plt.subplot( gs[ 1 ] )
# ax3 = plt.subplot( gs[ 2 ] )
# ax4 = plt.subplot( gs[ 3 ] )
# plotHandles = [ ax1, ax2, ax3, ax4 ]

trajectory_id_np = np.unique( trajectory_id )
# if len( trajectory_id_np ) > 20:
trajectory_set_1 = trajectory_id_np[ 0:-( len(trajectory_id_np)/2 ) ]
trajectory_set_2 = np.setdiff1d( trajectory_id_np, trajectory_set_1 )
color_set_1 = plt.cm.Vega10( np.linspace( 0, 1, len( trajectory_set_1 ) ) )
color_set_2 = plt.cm.Vega10( np.linspace( 0, 1, len( trajectory_set_2 ) ) )

## Plot altitude variation versus time now
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
for index in range( 0, len( trajectory_set_1 ) ):
    current_trajectory = trajectory_set_1[ index ]
    data_indices = np.where( trajectory_id == current_trajectory )
    data_indices = data_indices[ 0 ]

    xPlot = ( x[ data_indices ] )
    yPlot = ( y[ data_indices ] )
    zPlot = ( z[ data_indices ] )
    tPlot = ( t[ data_indices ] )
    tPlot = tPlot / ( 60.0 * 60.0 * 24.0 )
    rangePlot = np.sqrt( xPlot**2 + yPlot**2 + zPlot**2 )

    current_azimuth = launchAzimuth[ data_indices ].tolist( )
    current_azimuth = current_azimuth[ 0 ]
    current_velocity = velocityMagnitude[ data_indices ].tolist( )
    current_velocity = current_velocity[ 0 ]
    current_solarPhase = solarPhase[ data_indices ].tolist( )
    current_solarPhase = current_solarPhase[ 0 ]

    ax1.plot( tPlot, rangePlot, c=color_set_1[index],                                        \
              label='Launch azimuth = ' + str( current_azimuth ) + '[deg], '
                    + 'Launch Velocity = ' + str( current_velocity ) + '[m/s], '
                    + 'Solar phase = ' + str( current_solarPhase ) + '[deg]' )

ax1.axhline( y=2.0 * alpha, color='red', label='LAO' )
ax1.axhline( y=3.0 * alpha, color='blue', label='MAO' )
ax1.axhline( y=5.0 * alpha, color='orange', label='HAO' )
ax1.axhline( y=7.0 * alpha, color='green', label='UHAO' )

ax1.set_ylabel('Range [m]')
ax1.set_xlabel('Time [Days]')
ax1.set_yscale('log')
# ax1.set_xscale('log')
ax1.grid( True )
ax1.legend( ).draggable( )

fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax2 = plt.subplot( gs[ 0 ] )
for index in range( 0, len( trajectory_set_2 ) ):
    current_trajectory = trajectory_set_2[ index ]
    data_indices = np.where( trajectory_id == current_trajectory )
    data_indices = data_indices[ 0 ]

    xPlot = ( x[ data_indices ] )
    yPlot = ( y[ data_indices ] )
    zPlot = ( z[ data_indices ] )
    tPlot = ( t[ data_indices ] )
    tPlot = tPlot / ( 60.0 * 60.0 * 24.0 )
    rangePlot = np.sqrt( xPlot**2 + yPlot**2 + zPlot**2 )

    current_azimuth = launchAzimuth[ data_indices ].tolist( )
    current_azimuth = current_azimuth[ 0 ]
    current_velocity = velocityMagnitude[ data_indices ].tolist( )
    current_velocity = current_velocity[ 0 ]
    current_solarPhase = solarPhase[ data_indices ].tolist( )
    current_solarPhase = current_solarPhase[ 0 ]

    ax2.plot( tPlot, rangePlot, c=color_set_2[index],                                        \
              label='Launch azimuth = ' + str( current_azimuth ) + '[deg], '
                    + 'Launch Velocity = ' + str( current_velocity ) + '[m/s], '
                    + 'Solar phase = ' + str( current_solarPhase ) + '[deg]' )

ax2.axhline( y=2.0 * alpha, color='red', label='LAO' )
ax2.axhline( y=3.0 * alpha, color='blue', label='MAO' )
ax2.axhline( y=5.0 * alpha, color='orange', label='HAO' )
ax2.axhline( y=7.0 * alpha, color='green', label='UHAO' )

ax2.set_ylabel('Range [m]')
ax2.set_xlabel('Time [Days]')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.grid( True )
ax2.legend( ).draggable( )
# else:
#     fig = plt.figure( )
#     gs = gridspec.GridSpec( 1, 1 )
#     ax1 = plt.subplot( gs[ 0 ] )
#     colors = plt.cm.Vega20( np.linspace( 0, 1, len( trajectory_id_np ) ) )

#     ## Plot altitude variation versus time now
#     for index in range( 0, len( trajectory_id_np ) ):
#         current_trajectory = trajectory_id_np[ index ]
#         data_indices = np.where( trajectory_id == current_trajectory )
#         data_indices = data_indices[ 0 ]

#         xPlot = ( x[ data_indices ] )
#         yPlot = ( y[ data_indices ] )
#         zPlot = ( z[ data_indices ] )
#         tPlot = ( t[ data_indices ] )
#         tPlot = tPlot / ( 60.0 * 60.0 * 24.0 )
#         rangePlot = np.sqrt( xPlot**2 + yPlot**2 + zPlot**2 )

#         current_azimuth = launchAzimuth[ data_indices ].tolist( )
#         current_azimuth = current_azimuth[ 0 ]
#         current_velocity = velocityMagnitude[ data_indices ].tolist( )
#         current_velocity = current_velocity[ 0 ]
#         current_solarPhase = solarPhase[ data_indices ].tolist( )
#         current_solarPhase = current_solarPhase[ 0 ]

#         ax1.plot( tPlot, rangePlot, c=colors[index],                                        \
#                   label='Launch azimuth = ' + str( current_azimuth ) + '[deg], '
#                         + 'Launch Velocity = ' + str( current_velocity ) + '[m/s], '
#                         + 'Solar phase = ' + str( current_solarPhase ) + '[deg]' )

#     ax1.axhline( y=2.0 * alpha, color='red', label='LAO' )
#     ax1.axhline( y=3.0 * alpha, color='blue', label='MAO' )
#     ax1.axhline( y=5.0 * alpha, color='orange', label='HAO' )
#     ax1.axhline( y=7.0 * alpha, color='green', label='UHAO' )

#     ax1.set_ylabel('Altitude [m]')
#     ax1.set_xlabel('Time [Days]')
#     ax1.set_yscale('log')
#     ax1.set_xscale('log')
#     ax1.grid( True )
#     ax1.legend( ).draggable( )

# uniqueSolarPhase = np.unique( solarPhase )
# for index in range( 0, len( uniqueSolarPhase ) ):
#     print "Loop count = " + str( index + 1 )
#     currentSolarPhase = uniqueSolarPhase[ index ]
#     currentDataIndices = np.where( solarPhase == currentSolarPhase )
#     currentDataIndices = currentDataIndices[ 0 ]
#     currentPlotHandle = plotHandles[ index ]

#     velocities = velocityMagnitude[ currentDataIndices ]
#     azimuthsFirstExtract = launchAzimuth[ currentDataIndices ]
#     currentUniqueVelocities = np.unique( velocities )
#     for index2 in range( 0, len( currentUniqueVelocities ) ):
#         currentVelocity = currentUniqueVelocities[ index2 ]
#         currentDataIndices = np.where( velocities == currentVelocity )
#         currentDataIndices = currentDataIndices[ 0 ]

#         azimuthsSecondExtract = azimuthsFirstExtract[ currentDataIndices ]
#         currentUniqueAzimuths = np.unique( azimuthsSecondExtract )
#         for index3 in range( 0, len( currentUniqueAzimuths ) ):
#             currentAzimuth = currentUniqueAzimuths[ index3 ]
#             currentDataIndices = np.where( azimuthsSecondExtract == currentAzimuth )
#             currentDataIndices = currentDataIndices[ 0 ]

#             xPlot = x[ currentDataIndices ]
#             yPlot = y[ currentDataIndices ]
#             zPlot = z[ currentDataIndices ]
#             tPlot = t[ currentDataIndices ]
#             tPlot = tPlot / ( 60.0 * 60.0 * 24.0 )
#             rangePlot = np.sqrt( xPlot**2 + yPlot**2 + zPlot**2 )
#             currentPlotHandle.plot( tPlot, rangePlot,                     \
#                                     label='Launch azimuth = '
#                                     + str( currentAzimuth ) + '[deg], '
#                                     + 'Launch velocity = ' + str( currentVelocity ) + ' [m/s]' )

#     currentPlotHandle.set_yscale('log')
#     currentPlotHandle.set_xscale('log')
#     currentPlotHandle.set_xlabel('Time [days]')
#     currentPlotHandle.set_ylabel('Altitude [m]')
#     currentPlotHandle.axhline( y=2.0 * alpha, color='red', label='LAO' )
#     currentPlotHandle.axhline( y=3.0 * alpha, color='blue', label='MAO' )
#     currentPlotHandle.axhline( y=5.0 * alpha, color='orange', label='HAO' )
#     currentPlotHandle.axhline( y=7.0 * alpha, color='green', label='UHAO' )
#     currentPlotHandle.legend( ).draggable( )
#     currentPlotHandle.set_title( 'Solar phase angle = ' + str( currentSolarPhase ) + ' [deg]' )

## Unique azimuth
# unique_azimuths = np.unique( launchAzimuth )
# if len( unique_azimuths ) > 20:
#     unique_azimuths_set_1 = unique_azimuths[ 0:-( len(unique_azimuths)/2 ) ]
#     unique_azimuths_set_2 = np.setdiff1d( unique_azimuths, unique_azimuths_set_1 )
#     fig = plt.figure( )
#     gs = gridspec.GridSpec( 2, 1 )
#     ax1 = plt.subplot( gs[ 0 ] )
#     ax2 = plt.subplot( gs[ 1 ] )
#     color_set_1 = plt.cm.Vega20( np.linspace( 0, 1, len( unique_azimuths_set_1 ) ) )
#     color_set_2 = plt.cm.Vega20( np.linspace( 0, 1, len( unique_azimuths_set_2 ) ) )

#     ## Plot altitude variation versus time now
#     for index in range( 0, len( unique_azimuths_set_1 ) ):
#         current_azimuth = unique_azimuths_set_1[ index ]
#         data_indices = np.where( launchAzimuth == current_azimuth )
#         data_indices = data_indices[ 0 ]
#         xPlot = ( x[ data_indices ] )
#         yPlot = ( y[ data_indices ] )
#         zPlot = ( z[ data_indices ] )
#         tPlot = ( t[ data_indices ] )
#         tPlot = tPlot / ( 60.0 * 60.0 )
#         rangePlot = np.sqrt( xPlot**2 + yPlot**2 + zPlot**2 )
#         ax1.plot( tPlot, rangePlot, c=color_set_1[index],                                        \
#                   label='Launch azimuth = ' + str( current_azimuth ) + '[deg]' )

#     for index in range( 0, len( unique_azimuths_set_2 ) ):
#         current_azimuth = unique_azimuths_set_2[ index ]
#         data_indices = np.where( launchAzimuth == current_azimuth )
#         data_indices = data_indices[ 0 ]
#         xPlot = ( x[ data_indices ] )
#         yPlot = ( y[ data_indices ] )
#         zPlot = ( z[ data_indices ] )
#         tPlot = ( t[ data_indices ] )
#         tPlot = tPlot / ( 60.0 * 60.0 )
#         rangePlot = np.sqrt( xPlot**2 + yPlot**2 + zPlot**2 )
#         ax2.plot( tPlot, rangePlot, c=color_set_2[index],                                        \
#                   label='Launch azimuth = ' + str( current_azimuth ) + '[deg]' )

#     ax1.axhline( y=2.0 * alpha, color='red', label='LAO' )
#     ax1.axhline( y=3.0 * alpha, color='blue', label='MAO' )
#     ax1.axhline( y=5.0 * alpha, color='orange', label='HAO' )
#     ax1.axhline( y=7.0 * alpha, color='green', label='UHAO' )

#     ax1.set_ylabel('Altitude [m]')
#     ax1.set_xlabel('Time [hrs]')
#     ax1.set_yscale('log')
#     ax1.set_xscale('log')
#     # ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
#     ax1.grid( True )
#     ax1.legend( markerscale=7 ).draggable( )

#     ax2.axhline( y=2.0 * alpha, color='red', label='LAO' )
#     ax2.axhline( y=3.0 * alpha, color='blue', label='MAO' )
#     ax2.axhline( y=5.0 * alpha, color='orange', label='HAO' )
#     ax2.axhline( y=7.0 * alpha, color='green', label='UHAO' )

#     ax2.set_ylabel('Altitude [m]')
#     ax2.set_xlabel('Time [hrs]')
#     ax2.set_yscale('log')
#     ax2.set_xscale('log')
#     # ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
#     ax2.grid( True )
#     ax2.legend( markerscale=7 ).draggable( )
# else:
#     fig = plt.figure( )
#     gs = gridspec.GridSpec( 1, 1 )
#     ax1 = plt.subplot( gs[ 0 ] )
#     colors = plt.cm.Vega20( np.linspace( 0, 1, len( unique_azimuths ) ) )

#     ## Plot altitude variation versus time now
#     for index in range( 0, len( unique_azimuths ) ):
#         current_azimuth = unique_azimuths[ index ]
#         data_indices = np.where( launchAzimuth == current_azimuth )
#         data_indices = data_indices[ 0 ]
#         xPlot = ( x[ data_indices ] )
#         yPlot = ( y[ data_indices ] )
#         zPlot = ( z[ data_indices ] )
#         tPlot = ( t[ data_indices ] )
#         tPlot = tPlot / ( 60.0 * 60.0 )
#         rangePlot = np.sqrt( xPlot**2 + yPlot**2 + zPlot**2 )
#         ax1.plot( tPlot, rangePlot, c=colors[index],                                  \
#                   label='Launch azimuth = ' + str( current_azimuth ) + '[deg]' )

#     ax1.axhline( y=2.0 * alpha, color='red', label='LAO' )
#     ax1.axhline( y=3.0 * alpha, color='blue', label='MAO' )
#     ax1.axhline( y=5.0 * alpha, color='orange', label='HAO' )
#     ax1.axhline( y=7.0 * alpha, color='green', label='UHAO' )

#     ax1.set_ylabel('Altitude [m]')
#     ax1.set_xlabel('Time [hrs]')
#     ax1.set_yscale('log')
#     ax1.set_xscale('log')
#     # ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
#     ax1.grid( True )
#     ax1.legend( markerscale=7 ).draggable( )

# plt.suptitle( 'Time spent by Regolith in different altitude regimes \n'
#                + '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '
#                + 'Longest edge, Capture case' )

if database:
    database.close( )

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
