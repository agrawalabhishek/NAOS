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

maxTime = 9.0 * 30.0 * 24.0 * 60.0 * 60.0

phaseAngle = 'N.A.'

## Operations
# Connect to SQLite database.
try:
        # database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity_with_perturbations/phase_0/leadingEdge.db")
        # database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )
        database = sqlite3.connect( "../data/guarantee_escape_speed/longest_edge/normalLaunch.db" )
        # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

## get the local directional escape speeds in rotating and inertial frame for each launch azimuth
data4 = pd.read_sql( "SELECT    directional_escape_speed,                                   \
                                directional_inertial_escape_speed,                          \
                                ROUND( launch_declination ),                                \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( start_flag = 1 );",                                       \
                     database )

data4.columns = [ 'directional_escape_speed',                                               \
                  'directional_inertial_escape_speed',                                      \
                  'directional_escape_declination',                                         \
                  'directional_escape_azimuth' ]

directionalEscapeSpeed          = data4[ 'directional_escape_speed' ]
inertialDirectionalEscapeSpeed  = data4[ 'directional_inertial_escape_speed' ]
directionalEscapeDeclination    = data4[ 'directional_escape_declination' ]
directionalEscapeAzimuth        = data4[ 'directional_escape_azimuth' ]

## get data for escape cases
data1 = pd.read_sql( "SELECT    initial_velocity_x,                                         \
                                initial_velocity_y,                                         \
                                initial_velocity_z,                                         \
                                initial_inertial_velocity_x,                                \
                                initial_inertial_velocity_y,                                \
                                initial_inertial_velocity_z,                                \
                                ROUND( launch_declination ),                                \
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
                  'launch_declination',                                                     \
                  'launch_azimuth' ]

escape_rotFrame_vx                     = data1[ 'rotFrame_vx' ]
escape_rotFrame_vy                     = data1[ 'rotFrame_vy' ]
escape_rotFrame_vz                     = data1[ 'rotFrame_vz' ]
escape_inertial_vx                     = data1[ 'inertial_vx' ]
escape_inertial_vy                     = data1[ 'inertial_vy' ]
escape_inertial_vz                     = data1[ 'inertial_vz' ]
escape_declination                     = data1[ 'launch_declination' ]
escape_azimuth                         = data1[ 'launch_azimuth' ]

# Set up the figure
fig = plt.figure( )
plt.suptitle( 'Regolith launch and guaranteed escape velocity for a spherical asteroid\nEscape cases' )
gs = gridspec.GridSpec( 2, 1, height_ratios = [ 1, 1 ] )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )

## plot launch azimuth versus rotating frame initial velocity for escape cases
escape_rotFrameInitialVelocity = np.sqrt( escape_rotFrame_vx**2 + escape_rotFrame_vy**2 + escape_rotFrame_vz**2 )
ax1Handle1 = ax1.scatter( escape_azimuth, escape_rotFrameInitialVelocity, s=5, color=colors.cnames['purple'], label='actual speed' )
ax1Handle2, = ax1.plot( directionalEscapeAzimuth, directionalEscapeSpeed, color=colors.cnames['red'], label='Conservative escape speed', lw=2 )

## format axis and title
ax1.set_xlabel('Launch azimuth [deg]')
ax1.set_ylabel('$V_{initial}$ [m/s]')
ax1.set_xlim( 0, 360 )
# ax1.set_title( 'Regolith (escape) initial velocity versus launch azimuth (rotating frame) \n Phase angle = ' + str(phaseAngle) + '[deg]' )
ax1.set_title( 'Rotating frame' )
ax1.legend( ).draggable( )
ax1.grid( )

## plot launch azimuth versus inertial frame initial velocity for escape cases
escape_inertialInitialVelocity = np.sqrt( escape_inertial_vx**2 + escape_inertial_vy**2 + escape_inertial_vz**2 )
ax2Handle1 = ax2.scatter( escape_azimuth, escape_inertialInitialVelocity, s=5, color=colors.cnames['orange'], label='actual speed' )
ax2Handle2, = ax2.plot( directionalEscapeAzimuth, inertialDirectionalEscapeSpeed, color=colors.cnames['red'], label='Conservative escape speed', lw=2 )

## format axis and title
ax2.set_xlabel('Launch azimuth [deg]')
ax2.set_ylabel('$V_{initial}$ [m/s]')
ax2.set_xlim( 0, 360 )
ax2.set_title( 'Inertial frame' )
ax2.legend( ).draggable( )
ax2.grid( )

# plt.show( )
# sys.exit( )

## Plot launch azimuth versus inertial initial velocity for escape cases
# Set up the figure
fig = plt.figure( )
plt.suptitle( 'Inertial launch velocities for Escape and Re-impact cases' )
gs = gridspec.GridSpec( 2, 1, height_ratios = [ 1, 1 ] )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )

escape_inertialInitialVelocity = np.sqrt( escape_inertial_vx**2 + escape_inertial_vy**2 + escape_inertial_vz**2 )
# ax1Handle1 = ax1.scatter( escape_azimuth, escape_inertialInitialVelocity, s=5, color=colors.cnames['orange'], label='actual speed' )
ax1Handle1 = ax1.scatter( escape_declination, escape_inertialInitialVelocity, s=5, color=colors.cnames['orange'], label='actual speed' )
ax1Handle2, = ax1.plot( directionalEscapeDeclination, inertialDirectionalEscapeSpeed, color=colors.cnames['red'], label='non-conservative escape speed', lw=2 )

## format axis and title
ax1.set_xlabel('Launch declination [deg]')
ax1.set_ylabel('$V_{initial}$ [m/s]')
# ax1.set_xlim( 0, 360 )
ax1.set_title( 'Escape' )
ax1.legend( ).draggable( )
ax1.grid( )

## get data for re-impact cases
data2 = pd.read_sql( "SELECT    initial_velocity_x,                                         \
                                initial_velocity_y,                                         \
                                initial_velocity_z,                                         \
                                initial_inertial_velocity_x,                                \
                                initial_inertial_velocity_y,                                \
                                initial_inertial_velocity_z,                                \
                                ROUND( launch_declination ),                                \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( crash_flag = 1 );",                                       \
                     database )

data2.columns = [ 'rotFrame_vx',                                                            \
                  'rotFrame_vy',                                                            \
                  'rotFrame_vz',                                                            \
                  'inertial_vx',                                                            \
                  'inertial_vy',                                                            \
                  'inertial_vz',                                                            \
                  'launch_declination',                                                     \
                  'launch_azimuth' ]

crash_rotFrame_vx                     = data2[ 'rotFrame_vx' ]
crash_rotFrame_vy                     = data2[ 'rotFrame_vy' ]
crash_rotFrame_vz                     = data2[ 'rotFrame_vz' ]
crash_inertial_vx                     = data2[ 'inertial_vx' ]
crash_inertial_vy                     = data2[ 'inertial_vy' ]
crash_inertial_vz                     = data2[ 'inertial_vz' ]
crash_declination                     = data2[ 'launch_declination' ]
crash_azimuth                         = data2[ 'launch_azimuth' ]

## plot launch azimuth versus rot frame initial velocity for crash cases
# crash_rotFrameInitialVelocity = np.sqrt( crash_rotFrame_vx**2 + crash_rotFrame_vy**2 + crash_rotFrame_vz**2 )
# ax1Handle1 = ax1.scatter( crash_azimuth, crash_rotFrameInitialVelocity, color=colors.cnames['purple'], label='actual speed' )
# ax1Handle2, = ax1.plot( directionalEscapeAzimuth, directionalEscapeSpeed, color=colors.cnames['red'], label='conservative escape speed', lw=2 )

## format axis and title
# ax1.set_xlabel('Launch azimuth [deg]')
# ax1.set_ylabel('$V_{initial}$ [m/s]')
# ax1.set_xlim( 0, 360 )
# ax1.set_title( 'Regolith (re-impact) initial velocity versus launch azimuth (rotating frame) \n Phase angle = ' + str(phaseAngle) + '[deg]' )
# ax1.legend( ).draggable( )
# ax1.grid( )

## Plot launch azimuth versus inertial initial velocity for crash cases
crash_inertialInitialVelocity = np.sqrt( crash_inertial_vx**2 + crash_inertial_vy**2 + crash_inertial_vz**2 )
# ax2Handle1 = ax2.scatter( crash_azimuth, crash_inertialInitialVelocity, s=5, color=colors.cnames['purple'], label='actual speed' )
ax2Handle1 = ax2.scatter( crash_declination, crash_inertialInitialVelocity, s=5, color=colors.cnames['purple'], label='actual speed' )
ax2Handle2, = ax2.plot( directionalEscapeDeclination, inertialDirectionalEscapeSpeed, color=colors.cnames['red'], label='non-conservative escape speed', lw=2 )

## format axis and title
ax2.set_xlabel('Launch azimuth [deg]')
ax2.set_ylabel('$V_{initial}$ [m/s]')
# ax2.set_xlim( 0, 360 )
ax2.set_title( 'Reimpact' )
ax2.legend( ).draggable( )
ax2.grid( )

## get data for temporary capture cases
data3 = pd.read_sql( "SELECT    initial_velocity_x,                                         \
                                initial_velocity_y,                                         \
                                initial_velocity_z,                                         \
                                initial_inertial_velocity_x,                                \
                                initial_inertial_velocity_y,                                \
                                initial_inertial_velocity_z,                                \
                                ROUND( launch_azimuth ),                                    \
                                total_energy,                                               \
                                eccentricity                                                \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( end_flag = 1 );",                                         \
                     database )

data3.columns = [ 'rotFrame_vx',                                                            \
                  'rotFrame_vy',                                                            \
                  'rotFrame_vz',                                                            \
                  'inertial_vx',                                                            \
                  'inertial_vy',                                                            \
                  'inertial_vz',                                                            \
                  'launch_azimuth',                                                         \
                  'energy',                                                                 \
                  'eccentricity' ]

capture_rotFrame_vx                     = data3[ 'rotFrame_vx' ]
capture_rotFrame_vy                     = data3[ 'rotFrame_vy' ]
capture_rotFrame_vz                     = data3[ 'rotFrame_vz' ]
capture_inertial_vx                     = data3[ 'inertial_vx' ]
capture_inertial_vy                     = data3[ 'inertial_vy' ]
capture_inertial_vz                     = data3[ 'inertial_vz' ]
capture_azimuth                         = data3[ 'launch_azimuth' ]
capture_energy                          = data3[ 'energy' ]
capture_eccentricity                    = data3[ 'eccentricity' ]

# limit eccentricity decimal digits for a cleaner display
for i in range( 0, (len(capture_eccentricity)) ):
    capture_eccentricity[ i ] = float( "{0:.5f}".format( capture_eccentricity[ i ] ) )

capture_inertialInitialVelocity = np.sqrt( capture_inertial_vx**2 + capture_inertial_vy**2 + capture_inertial_vz**2 )
capture_rotFrameInitialVelocity = np.sqrt( capture_rotFrame_vx**2 + capture_rotFrame_vy**2 + capture_rotFrame_vz**2 )

if capture_azimuth.size != 0:
    ## Set up the figure
    # fig = plt.figure( )
    # gs = gridspec.GridSpec( 3, 1, height_ratios = [ 1, 1, 2.5 ] )
    # ax1 = plt.subplot( gs[ 0 ] )
    # ax2 = plt.subplot( gs[ 1 ] )
    # ax3 = plt.subplot( gs[ 2 ], frameon=False )

    # ## plot launch azimuth versus rot frame initial velocity for capture cases
    # ax1Handle1 = ax1.scatter( capture_azimuth, capture_rotFrameInitialVelocity,
    #                           color=colors.cnames['purple'],
    #                           label='actual speed' )

    # ## add text box indicating energy and eccentricity at the scatter points
    # # for i, txt in enumerate( capture_eccentricity ):
    # #     ax1.annotate( "e="+str(txt),                                                                    \
    # #                   xy=( capture_azimuth[ i ], capture_rotFrameInitialVelocity[ i ] ),                \
    # #                   xytext=( capture_azimuth[ i ] - 10, capture_rotFrameInitialVelocity[ i ] + 1 ),   \
    # #                   arrowprops=dict( facecolor='black', shrink=0.05 ) )

    # ## Plot launch azimuth versus inertial initial velocity for capture cases
    # ax2Handle1 = ax2.scatter( capture_azimuth, capture_inertialInitialVelocity, color=colors.cnames['purple'], label='actual speed' )

    # # ax1Handle2, = ax1.plot( directionalEscapeAzimuth, directionalEscapeSpeed, color=colors.cnames['red'], label='directional escape speed', lw=2 )
    # # ax2Handle2, = ax2.plot( directionalEscapeAzimuth, inertialDirectionalEscapeSpeed, color=colors.cnames['red'], label='directional escape speed', lw=2 )

    # ## format axis and title
    # ax1.set_xlabel('Launch azimuth [deg]')
    # ax1.set_ylabel('$V_{initial}$ [m/s]')
    # ax1.set_xlim( 0, 360 )
    # ax1.set_title( 'Regolith (capture) initial velocity versus launch azimuth (rotating frame) \n Ellipsoidal asteroid, Phase angle = ' + str(phaseAngle) + '[deg]' )
    # ax1.legend( )
    # # ax1.legend( loc='center left', bbox_to_anchor=( 0.995,0.5 ) )
    # ax1.grid( )

    # ## format axis and title
    # ax2.set_xlabel('Launch azimuth [deg]')
    # ax2.set_ylabel('$V_{initial}$ [m/s]')
    # ax2.set_xlim( 0, 360 )
    # ax2.set_title( 'Regolith (capture) initial velocity versus launch azimuth (inertial frame) \n Ellipsoidal asteroid, Phase angle = ' + str(phaseAngle) + '[deg]' )
    # ax2.legend( )
    # # ax2.legend( loc='center left', bbox_to_anchor=( 0.995,0.5 ) )
    # ax2.grid( )

    fig = plt.figure( )
    gs = gridspec.GridSpec( 1, 1 )
    ax3 = plt.subplot( gs[ 0 ], frameon=False )

    ## add the metatable
    ax3.axis( 'off' )
    metadata_table = []
    # metadata_table.append( [ "velocity (rotating Frame) [m/s]", \
    #                          "launch azimuth [deg]",            \
    #                          "energy [m^2 / s^2]",              \
    #                          "eccentricity" ] )
    columnLabels = [ "V (rotating Frame) [m/s]",        \
                     "launch azimuth [deg]",            \
                     "energy [$m^2 / s^2$]",            \
                     "eccentricity" ]

    for i in range( 0, len( capture_energy ) ):
        metadata_table.append( [ capture_rotFrameInitialVelocity[ i ],  \
                                 capture_azimuth[ i ],                  \
                                 capture_energy[ i ],                   \
                                 capture_eccentricity[ i ] ] )

    table = ax3.table( cellText = metadata_table, colLabels = columnLabels, cellLoc = 'center', loc = 'center' )
    table.auto_set_font_size(False)
    table.set_fontsize( 14 )
    table_properties = table.properties( )
    table_cells = table_properties[ 'child_artists' ]
    for cell in table_cells: cell.set_height( 0.15 )
    cell_dict = table.get_celld( )
    # for row in xrange( 0, len( capture_energy ) ): cell_dict[ ( row, 2 ) ].set_width( 0.1 )

if database:
    database.close( )

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
