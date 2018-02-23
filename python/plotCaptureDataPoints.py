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

# Connect to SQLite database for the perturbing case
try:
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "3.2Density_1cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data for the perturbing case now...\n"

data1 = pd.read_sql( "SELECT    trajectory_id                                       \
                      FROM      regolith_trajectory_results                         \
                      WHERE     end_flag = 1;",                                     \
                      database )

data1.columns = [ 'trajectory_id' ]

trajectory_id = data1[ 'trajectory_id' ]
trajectory_id = tuple( trajectory_id.tolist( ) )

data2 = pd.read_sql( "SELECT    ROUND( initial_solar_phase_angle ),                     \
                                ROUND( initial_velocity_magnitude ),                    \
                                ROUND( launch_azimuth )                                 \
                      FROM      regolith_trajectory_results                             \
                      WHERE     start_flag = 1                                          \
                      AND       trajectory_id IN " + str( trajectory_id ) + ";",        \
                      database )

data2.columns = [ 'solar_phase_angle',                                                  \
                  'launch_velocity',                                                    \
                  'launch_azimuth' ]

solarPhase      = data2[ 'solar_phase_angle' ]
launchVelocity  = data2[ 'launch_velocity' ]
launchAzimuth   = data2[ 'launch_azimuth' ]

extractVelocities   = tuple( launchVelocity.tolist( ) )
extractAzimuths     = tuple( launchAzimuth.tolist( ) )

print "\nCapture velocities = " + str( extractVelocities )
print "\nCapture azimuths = " + str( extractAzimuths )
print "\nTotal capture cases = " + str( len( extractVelocities ) )

if database:
    database.close( )

try:
    database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                               + "multiple_launch_velocity/"
                               + "simulation_time_9_months/"
                               + "longestEdge.db")

except sqlite3.Error, e:
    print "Error %s:" % e.args[0]
    sys.exit(1)

print "Extracting data for the non-perturbing case now...\n"

data3 = pd.read_sql( "SELECT    trajectory_id,                                                  \
                                escape_flag,                                                    \
                                crash_flag,                                                     \
                                end_flag,                                                       \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth )                                         \
                      FROM      regolith_trajectory_results                                     \
                      WHERE     ROUND( initial_velocity_magnitude ) IN "
                                + str( extractVelocities ) + " " +
                     "AND       ROUND( launch_azimuth ) IN "
                                + str( extractAzimuths ) + " ;",
                     database )

data3.columns = [ 'trajectory_id',                                                              \
                  'escapeFlag',                                                                 \
                  'crashFlag',                                                                  \
                  'endFlag',                                                                    \
                  'launch_velocity',                                                            \
                  'launch_azimuth' ]

nonPerturbingTrajectoryId   = data3[ 'trajectory_id' ]
escapeFlag                  = data3[ 'escapeFlag' ]
crashFlag                   = data3[ 'crashFlag' ]
endFlag                     = data3[ 'endFlag' ]
nonPerturbingVelocity       = data3[ 'launch_velocity' ]
nonPerturbingAzimuth        = data3[ 'launch_azimuth' ]

print "Processing all data now...\n"

## extract information on the final state of regolith in non perturbing simulations for the
## corresponding capture cases in the perturbing simulation
uniqueNonPerturbingVelocities = np.unique( nonPerturbingVelocity )
uniqueNonPerturbingTrajectories = np.unique( nonPerturbingTrajectoryId )

nonPerturbingEscapeCase     = np.zeros( ( len( extractVelocities ), 2 ) )
nonPerturbingCrashCase      = np.zeros( ( len( extractVelocities ), 2 ) )
nonPerturbingCaptureCase    = np.zeros( ( len( extractVelocities ), 2 ) )
storingCounterEscape = 0
storingCounterCrash = 0
storingCounterCapture = 0

for index in range( 0, len( uniqueNonPerturbingTrajectories ) ):
    currentTrajectory           = uniqueNonPerturbingTrajectories[ index ]
    currentTrajectoryIndices    = np.where( nonPerturbingTrajectoryId == currentTrajectory )
    currentTrajectoryIndices    = currentTrajectoryIndices[ 0 ]
    currentDataSetLength        = len( currentTrajectoryIndices )

    currentLaunchVelocity = nonPerturbingVelocity[ currentTrajectoryIndices ].tolist( )
    currentLaunchVelocity = currentLaunchVelocity[ currentDataSetLength - 1 ]

    currentLaunchAzimuth = nonPerturbingAzimuth[ currentTrajectoryIndices ].tolist( )
    currentLaunchAzimuth = currentLaunchAzimuth[ currentDataSetLength - 1 ]

    for secondIndex in range( 0, len( launchVelocity ) ):
        if currentLaunchVelocity == launchVelocity[ secondIndex ] and currentLaunchAzimuth == launchAzimuth[ secondIndex ]:
            currentEscapeFlagValues = escapeFlag[ currentTrajectoryIndices ].tolist( )
            currentEscapeFlagValue = currentEscapeFlagValues[ currentDataSetLength - 1 ]
            if currentEscapeFlagValue == 1:
                nonPerturbingEscapeCase[ storingCounterEscape ][ 0 ] = currentLaunchVelocity
                nonPerturbingEscapeCase[ storingCounterEscape ][ 1 ] = currentLaunchAzimuth
                storingCounterEscape = storingCounterEscape + 1
                break

            currentCrashFlagValues = crashFlag[ currentTrajectoryIndices ].tolist( )
            currentCrashFlagValue = currentCrashFlagValues[ currentDataSetLength - 1 ]
            if currentCrashFlagValue == 1:
                nonPerturbingCrashCase[ storingCounterCrash ][ 0 ] = currentLaunchVelocity
                nonPerturbingCrashCase[ storingCounterCrash ][ 1 ] = currentLaunchAzimuth
                storingCounterCrash = storingCounterCrash + 1
                break

            currentEndFlagValues = endFlag[ currentTrajectoryIndices ].tolist( )
            currentEndFlagValue = currentEndFlagValues[ currentDataSetLength - 1 ]
            if currentEndFlagValue == 1:
                nonPerturbingCaptureCase[ storingCounterCapture ][ 0 ] = currentLaunchVelocity
                nonPerturbingCaptureCase[ storingCounterCapture ][ 1 ] = currentLaunchAzimuth
                storingCounterCapture = storingCounterCapture + 1
                break

## Set up the figure
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

## metadata table
ax1.axis( 'off' )
metadata_table = []
columnLabels = [ "Launch azimuth [deg]",                            \
                 "Launch velocity [m/s]",                           \
                 "Solar phase angle [deg]",                         \
                 "Non-perturbing status" ]

for i in range( 0, len( solarPhase ) ):
    currentAzimuth = launchAzimuth[ i ]
    currentVelocity = launchVelocity[ i ]

    for index in range( 0, len( nonPerturbingEscapeCase ) ):
        if currentAzimuth == nonPerturbingEscapeCase[ index ][ 1 ] and currentVelocity == nonPerturbingEscapeCase[ index ][ 0 ]:
            metadata_table.append( [ launchAzimuth[ i ],                    \
                                     launchVelocity[ i ],                   \
                                     solarPhase[ i ],                       \
                                     "Escape" ] )
            continue

    for index in range( 0, len( nonPerturbingCrashCase ) ):
        if currentAzimuth == nonPerturbingCrashCase[ index ][ 1 ] and currentVelocity == nonPerturbingCrashCase[ index ][ 0 ]:
            metadata_table.append( [ launchAzimuth[ i ],                    \
                                     launchVelocity[ i ],                   \
                                     solarPhase[ i ],                       \
                                     "Re-impact" ] )
            continue

    for index in range( 0, len( nonPerturbingCaptureCase ) ):
        if currentAzimuth == nonPerturbingCaptureCase[ index ][ 1 ] and currentVelocity == nonPerturbingCaptureCase[ index ][ 0 ]:
            metadata_table.append( [ launchAzimuth[ i ],                    \
                                     launchVelocity[ i ],                   \
                                     solarPhase[ i ],                       \
                                     "Temporary capture" ] )
            continue

table = ax1.table( cellText = metadata_table,
                   colLabels = columnLabels,
                   cellLoc = 'center',
                   loc = 'center' )
table.auto_set_font_size(False)
table.set_fontsize( 10 )
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

plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""

