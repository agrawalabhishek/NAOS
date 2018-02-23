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

def extractDataFromSQL( databaseConnect ):
    # general function for data extraction
    # get data for escape cases
    data1 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude ),                        \
                                    ROUND( initial_solar_phase_angle ),                         \
                                    ROUND( launch_azimuth )                                     \
                         FROM       regolith_trajectory_results                                 \
                         WHERE      ( escape_flag = 1 );",                                      \
                         databaseConnect )

    data1.columns = [ 'vel_mag',                                                                \
                      'solar_phase_angle',                                                      \
                      'launch_azimuth' ]

    escape_velocity_magnitude              = data1[ 'vel_mag' ]
    escape_azimuth                         = data1[ 'launch_azimuth' ]
    escape_solarPhase                      = data1[ 'solar_phase_angle' ]

    ## get data for re-impact cases
    data2 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude ),                        \
                                    ROUND( initial_solar_phase_angle ),                         \
                                    ROUND( launch_azimuth )                                     \
                         FROM       regolith_trajectory_results                                 \
                         WHERE      ( crash_flag = 1 );",                                       \
                         databaseConnect )

    data2.columns = [ 'vel_mag',                                                                \
                      'solar_phase_angle',                                                      \
                      'launch_azimuth' ]

    crash_velocity_magnitude              = data2[ 'vel_mag' ]
    crash_azimuth                         = data2[ 'launch_azimuth' ]
    crash_solarPhase                      = data2[ 'solar_phase_angle' ]

    ## get data for temporary capture cases
    data3 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude ),                        \
                                    ROUND( initial_solar_phase_angle ),                         \
                                    ROUND( launch_azimuth )                                     \
                         FROM       regolith_trajectory_results                                 \
                         WHERE      ( end_flag = 1 );",                                         \
                         databaseConnect )

    data3.columns = [ 'vel_mag',                                                                \
                      'solar_phase_angle',                                                      \
                      'launch_azimuth' ]

    capture_velocity_magnitude              = data3[ 'vel_mag' ]
    capture_azimuth                         = data3[ 'launch_azimuth' ]
    capture_solarPhase                      = data3[ 'solar_phase_angle' ]

    if databaseConnect:
        databaseConnect.close( )

    return ( escape_velocity_magnitude, escape_azimuth, escape_solarPhase,
             crash_velocity_magnitude, crash_azimuth, crash_solarPhase,
             capture_velocity_magnitude, capture_azimuth, capture_solarPhase )

def areaToMassRatio( radius, density ):
    crossSectionalArea = np.pi * radius**2
    mass = ( 4.0 / 3.0 ) * np.pi * radius**3 * density
    ratio = crossSectionalArea / mass
    return ratio

# Start timer.
start_time = time.time( )

## Operations
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593
mu = 876514

uniqueSolarPhase = np.array([ 45.0, 135.0, 225.0, 315.0 ])

## Operations
# Connect to SQLite database.
print "Extracting data now...\n"

commonPath = "/media/abhishek/Ashish/Thesis_Simulation_Databases/trailing_edge_perturbations_CDE/"

try:
    database1 = sqlite3.connect(commonPath + "3.2Density_1cmRadius/trailingEdge_3P2Density_1cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case1_launchVel_escape, case1_azimuth_escape, case1_solarPhase_escape,
  case1_launchVel_reimpact, case1_azimuth_reimpact, case1_solarPhase_reimpact,
  case1_launchVel_capture, case1_azimuth_capture, case1_solarPhase_capture )                    \
                                                        = extractDataFromSQL( database1 )
ratio = areaToMassRatio( 1.0e-2, 3.2e3 )
label1 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

try:
    database2 = sqlite3.connect(commonPath + "7.5Density_1cmRadius/trailingEdge_7P5Density_1cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case2_launchVel_escape, case2_azimuth_escape, case2_solarPhase_escape,
  case2_launchVel_reimpact, case2_azimuth_reimpact, case2_solarPhase_reimpact,
  case2_launchVel_capture, case2_azimuth_capture, case2_solarPhase_capture )                    \
                                                        = extractDataFromSQL( database2 )

ratio = areaToMassRatio( 1.0e-2, 7.5e3 )
label2 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

try:
    database3 = sqlite3.connect(commonPath + "3.2Density_5cmRadius/trailingEdge_3P2Density_5cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case3_launchVel_escape, case3_azimuth_escape, case3_solarPhase_escape,
  case3_launchVel_reimpact, case3_azimuth_reimpact, case3_solarPhase_reimpact,
  case3_launchVel_capture, case3_azimuth_capture, case3_solarPhase_capture )                    \
                                                        = extractDataFromSQL( database3 )

ratio = areaToMassRatio( 5.0e-2, 3.2e3 )
label3 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

try:
    database4 = sqlite3.connect(commonPath + "7.5Density_5cmRadius/trailingEdge_7P5Density_5cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case4_launchVel_escape, case4_azimuth_escape, case4_solarPhase_escape,
  case4_launchVel_reimpact, case4_azimuth_reimpact, case4_solarPhase_reimpact,
  case4_launchVel_capture, case4_azimuth_capture, case4_solarPhase_capture )                    \
                                                        = extractDataFromSQL( database4 )

ratio = areaToMassRatio( 5.0e-2, 7.5e3 )
label4 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

# try:
#     database5 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
#                                 + "multiple_launch_velocity_with_perturbations/"
#                                 + "simulation_time_9_months/"
#                                 + "3.2Density_10cmRadius/longestEdgePerturbations.db")

# except sqlite3.Error, e:
#         print "Error %s:" % e.args[0]
#         sys.exit(1)

# ( case5_launchVel_escape, case5_azimuth_escape, case5_solarPhase_escape,
#   case5_launchVel_reimpact, case5_azimuth_reimpact, case5_solarPhase_reimpact,
#   case5_launchVel_capture, case5_azimuth_capture, case5_solarPhase_capture )                    \
#                                                         = extractDataFromSQL( database5 )

# ratio = areaToMassRatio( 10.0e-2, 3.2e3 )
# label5 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

# try:
#     database6 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
#                                 + "multiple_launch_velocity_with_perturbations/"
#                                 + "simulation_time_9_months/"
#                                 + "7.5Density_10cmRadius/longestEdgePerturbations.db")

# except sqlite3.Error, e:
#         print "Error %s:" % e.args[0]
#         sys.exit(1)

# ( case6_launchVel_escape, case6_azimuth_escape, case6_solarPhase_escape,
#   case6_launchVel_reimpact, case6_azimuth_reimpact, case6_solarPhase_reimpact,
#   case6_launchVel_capture, case6_azimuth_capture, case6_solarPhase_capture )                    \
#                                                         = extractDataFromSQL( database6 )

# ratio = areaToMassRatio( 10.0e-2, 7.5e3 )
# label6 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

all_distinct_velocities = 16
# barColors = ['orange', 'crimson', 'green', 'cyan', 'grey', 'purple']
barColors = ['orange', 'crimson', 'green', 'cyan']

####################################################################################################
# plot the histogram for number of particles ESCAPED against initial launch velocity
fig = plt.figure( figsize=(12, 12) )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )
plotHandles = [ ax1, ax2, ax3, ax4 ]

for index in range( 0, len( uniqueSolarPhase ) ):
    currentPlotHandle = plotHandles[ index ]
    currentSolarPhase = uniqueSolarPhase[ index ]

    case1EscapeDataIndices = np.where( case1_solarPhase_escape == currentSolarPhase )
    case1EscapeDataIndices = case1EscapeDataIndices[ 0 ]
    case1_launchVel_escape_plot = case1_launchVel_escape[ case1EscapeDataIndices ]

    case2EscapeDataIndices = np.where( case2_solarPhase_escape == currentSolarPhase )
    case2EscapeDataIndices = case2EscapeDataIndices[ 0 ]
    case2_launchVel_escape_plot = case2_launchVel_escape[ case2EscapeDataIndices ]

    case3EscapeDataIndices = np.where( case3_solarPhase_escape == currentSolarPhase )
    case3EscapeDataIndices = case3EscapeDataIndices[ 0 ]
    case3_launchVel_escape_plot = case3_launchVel_escape[ case3EscapeDataIndices ]

    case4EscapeDataIndices = np.where( case4_solarPhase_escape == currentSolarPhase )
    case4EscapeDataIndices = case4EscapeDataIndices[ 0 ]
    case4_launchVel_escape_plot = case4_launchVel_escape[ case4EscapeDataIndices ]

    # case5EscapeDataIndices = np.where( case5_solarPhase_escape == currentSolarPhase )
    # case5EscapeDataIndices = case5EscapeDataIndices[ 0 ]
    # case5_launchVel_escape_plot = case5_launchVel_escape[ case5EscapeDataIndices ]

    # case6EscapeDataIndices = np.where( case6_solarPhase_escape == currentSolarPhase )
    # case6EscapeDataIndices = case6EscapeDataIndices[ 0 ]
    # case6_launchVel_escape_plot = case6_launchVel_escape[ case6EscapeDataIndices ]

    currentPlotHandle.hist( [case1_launchVel_escape_plot, case2_launchVel_escape_plot,
                             case3_launchVel_escape_plot, case4_launchVel_escape_plot],
                            bins=range(1, all_distinct_velocities+2),
                            histtype='bar',
                            # stacked=True,
                            color=barColors,
                            label=[label1, label2, label3, label4] if index == 1 else "__nolegend__",
                            align='left', alpha=1.0 )

    currentPlotHandle.set_xlim( 0, all_distinct_velocities+1 )
    currentPlotHandle.xaxis.set_ticks( np.arange(0, all_distinct_velocities+1, 1) )
    currentPlotHandle.set_xlabel( '$V_{initial}$ [m/s]' )
    currentPlotHandle.set_ylabel( 'Number of particles' )
    currentPlotHandle.set_title( 'Solar phase angle = ' + str(currentSolarPhase) + '[deg]' )
    currentPlotHandle.grid(True)
    if index == 1:
        currentPlotHandle.legend( ).draggable( )

plt.suptitle( 'Escape case histogram, Ellipsoid trailing edge' )

####################################################################################################
# plot the histogram for number of particles REIMPACTED against initial launch velocity
fig = plt.figure( figsize=(12, 12) )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )
plotHandles = [ ax1, ax2, ax3, ax4 ]

for index in range( 0, len( uniqueSolarPhase ) ):
    currentPlotHandle = plotHandles[ index ]
    currentSolarPhase = uniqueSolarPhase[ index ]

    case1reimpactDataIndices = np.where( case1_solarPhase_reimpact == currentSolarPhase )
    case1reimpactDataIndices = case1reimpactDataIndices[ 0 ]
    case1_launchVel_reimpact_plot = case1_launchVel_reimpact[ case1reimpactDataIndices ]

    case2reimpactDataIndices = np.where( case2_solarPhase_reimpact == currentSolarPhase )
    case2reimpactDataIndices = case2reimpactDataIndices[ 0 ]
    case2_launchVel_reimpact_plot = case2_launchVel_reimpact[ case2reimpactDataIndices ]

    case3reimpactDataIndices = np.where( case3_solarPhase_reimpact == currentSolarPhase )
    case3reimpactDataIndices = case3reimpactDataIndices[ 0 ]
    case3_launchVel_reimpact_plot = case3_launchVel_reimpact[ case3reimpactDataIndices ]

    case4reimpactDataIndices = np.where( case4_solarPhase_reimpact == currentSolarPhase )
    case4reimpactDataIndices = case4reimpactDataIndices[ 0 ]
    case4_launchVel_reimpact_plot = case4_launchVel_reimpact[ case4reimpactDataIndices ]

    # case5reimpactDataIndices = np.where( case5_solarPhase_reimpact == currentSolarPhase )
    # case5reimpactDataIndices = case5reimpactDataIndices[ 0 ]
    # case5_launchVel_reimpact_plot = case5_launchVel_reimpact[ case5reimpactDataIndices ]

    # case6reimpactDataIndices = np.where( case6_solarPhase_reimpact == currentSolarPhase )
    # case6reimpactDataIndices = case6reimpactDataIndices[ 0 ]
    # case6_launchVel_reimpact_plot = case6_launchVel_reimpact[ case6reimpactDataIndices ]

    currentPlotHandle.hist( [case1_launchVel_reimpact_plot, case2_launchVel_reimpact_plot,
                             case3_launchVel_reimpact_plot, case4_launchVel_reimpact_plot],
                            bins=range(1, all_distinct_velocities+2),
                            histtype='bar',
                            # stacked=True,
                            color=barColors,
                            label=[label1, label2, label3, label4] if index == 1 else "__nolegend__",
                            align='left', alpha=1.0 )

    currentPlotHandle.set_xlim( 0, all_distinct_velocities+1 )
    currentPlotHandle.xaxis.set_ticks( np.arange(0, all_distinct_velocities+1, 1) )
    currentPlotHandle.set_xlabel( '$V_{initial}$ [m/s]' )
    currentPlotHandle.set_ylabel( 'Number of particles' )
    currentPlotHandle.set_title( 'Solar phase angle = ' + str(currentSolarPhase) + '[deg]' )
    currentPlotHandle.grid(True)
    if index == 1:
        currentPlotHandle.legend( ).draggable( )

plt.suptitle( 'Reimpact case histogram, Ellipsoid trailing edge' )

####################################################################################################
# plot the histogram for number of particles CAPTURED against initial launch velocity
fig = plt.figure( figsize=(12, 12) )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )
plotHandles = [ ax1, ax2, ax3, ax4 ]

for index in range( 0, len( uniqueSolarPhase ) ):
    currentPlotHandle = plotHandles[ index ]
    currentSolarPhase = uniqueSolarPhase[ index ]

    case1captureDataIndices = np.where( case1_solarPhase_capture == currentSolarPhase )
    case1captureDataIndices = case1captureDataIndices[ 0 ]
    case1_launchVel_capture_plot = case1_launchVel_capture[ case1captureDataIndices ]

    case2captureDataIndices = np.where( case2_solarPhase_capture == currentSolarPhase )
    case2captureDataIndices = case2captureDataIndices[ 0 ]
    case2_launchVel_capture_plot = case2_launchVel_capture[ case2captureDataIndices ]

    case3captureDataIndices = np.where( case3_solarPhase_capture == currentSolarPhase )
    case3captureDataIndices = case3captureDataIndices[ 0 ]
    case3_launchVel_capture_plot = case3_launchVel_capture[ case3captureDataIndices ]

    case4captureDataIndices = np.where( case4_solarPhase_capture == currentSolarPhase )
    case4captureDataIndices = case4captureDataIndices[ 0 ]
    case4_launchVel_capture_plot = case4_launchVel_capture[ case4captureDataIndices ]

    # case5captureDataIndices = np.where( case5_solarPhase_capture == currentSolarPhase )
    # case5captureDataIndices = case5captureDataIndices[ 0 ]
    # case5_launchVel_capture_plot = case5_launchVel_capture[ case5captureDataIndices ]

    # case6captureDataIndices = np.where( case6_solarPhase_capture == currentSolarPhase )
    # case6captureDataIndices = case6captureDataIndices[ 0 ]
    # case6_launchVel_capture_plot = case6_launchVel_capture[ case6captureDataIndices ]

    currentPlotHandle.hist( [case1_launchVel_capture_plot,
                             case2_launchVel_capture_plot,
                             case3_launchVel_capture_plot,
                             case4_launchVel_capture_plot],
                            bins=range(1, all_distinct_velocities+2),
                            histtype='bar',
                            stacked=True,
                            color=barColors,
                            label=[label1, label2, label3, label4] if index == 1 else "__nolegend__",
                            align='left', alpha=1.0 )

    currentPlotHandle.set_xlim( 0, all_distinct_velocities+1 )
    currentPlotHandle.xaxis.set_ticks( np.arange(0, all_distinct_velocities+1, 1) )
    currentPlotHandle.set_xlabel( '$V_{initial}$ [m/s]' )
    currentPlotHandle.set_ylabel( 'Number of particles' )
    currentPlotHandle.set_title( 'Solar phase angle = ' + str(currentSolarPhase) + '[deg]' )
    currentPlotHandle.grid(True)
    if index == 1:
        currentPlotHandle.legend( ).draggable( )

plt.suptitle( 'Capture case histogram, Ellipsoid trailing edge' )

####################################################################################################

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
