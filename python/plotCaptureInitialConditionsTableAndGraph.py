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
    data2 = pd.read_sql( "SELECT    ROUND( initial_solar_phase_angle ),                     \
                                    ROUND( initial_velocity_magnitude ),                    \
                                    ROUND( launch_azimuth )                                 \
                          FROM      regolith_trajectory_results                             \
                          WHERE     end_flag = 1;",                                         \
                          databaseConnect )

    data2.columns = [ 'solar_phase_angle',                                                  \
                      'launch_velocity',                                                    \
                      'launch_azimuth' ]

    solarPhase      = data2[ 'solar_phase_angle' ]
    launchVelocity  = data2[ 'launch_velocity' ]
    launchAzimuth   = data2[ 'launch_azimuth' ]

    if databaseConnect:
        databaseConnect.close( )

    return ( solarPhase, launchVelocity, launchAzimuth )

def areaToMassRatio( radius, density ):
    crossSectionalArea = np.pi * radius**2
    mass = ( 4.0 / 3.0 ) * np.pi * radius**3 * density
    ratio = crossSectionalArea / mass
    return ratio

def plotDataAndTable( solarPhase, launchVelocity, launchAzimuth ):
    ## Set up the figure for the table to display initial conditions for the capture cases
    fig = plt.figure( )
    gs = gridspec.GridSpec( 1, 1 )
    ax1 = plt.subplot( gs[ 0 ] )

    ## metadata table
    ax1.axis( 'off' )
    metadata_table = []
    columnLabels = [ "Launch azimuth [deg]",                            \
                     "Launch velocity [m/s]",                           \
                     "Solar phase angle [deg]" ]

    for i in range( 0, len( solarPhase ) ):
        metadata_table.append( [ launchAzimuth[ i ],                    \
                                 launchVelocity[ i ],                   \
                                 solarPhase[ i ] ] )

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

def singlePhasePlotDataExtraction( currentSolarPhase,
                                   case_solarPhase, case_launchVelocity, case_launchAzimuth ):
    currentDataIndices = np.where( case_solarPhase == currentSolarPhase )
    currentDataIndices = currentDataIndices[ 0 ]
    currentLaunchVelocities = case_launchVelocity[ currentDataIndices ]
    currentLaunchAzimuths = case_launchAzimuth[ currentDataIndices ]
    return ( currentLaunchVelocities, currentLaunchAzimuths )

def plotCompartiveCaptureCases( case1_solarPhase, case1_launchVelocity, case1_launchAzimuth, label1,
                                case2_solarPhase, case2_launchVelocity, case2_launchAzimuth, label2,
                                case3_solarPhase, case3_launchVelocity, case3_launchAzimuth, label3,
                                case4_solarPhase, case4_launchVelocity, case4_launchAzimuth, label4,
                                case5_solarPhase, case5_launchVelocity, case5_launchAzimuth, label5,
                                case6_solarPhase, case6_launchVelocity, case6_launchAzimuth, label6 ):
    ## set up the figure to plot the data (solar phase wise)
    fig = plt.figure( figsize=(12, 12) )
    gs = gridspec.GridSpec( 2, 2 )
    ax1 = plt.subplot( gs[ 0 ] )
    ax2 = plt.subplot( gs[ 1 ] )
    ax3 = plt.subplot( gs[ 2 ] )
    ax4 = plt.subplot( gs[ 3 ] )

    plotHandles = [ ax1, ax2, ax3, ax4 ]
    uniqueSolarPhases = np.arange( 45.0, 360.0, 90.0 )
    for index in range( 0, len( uniqueSolarPhases ) ):
        currentPlotHandle = plotHandles[ index ]
        currentSolarPhase = uniqueSolarPhases[ index ]

        ( case1_currentLaunchVelocities, case1_currentLaunchAzimuth )                           \
                = singlePhasePlotDataExtraction( currentSolarPhase,
                                                 case1_solarPhase,
                                                 case1_launchVelocity,
                                                 case1_launchAzimuth )

        ( case2_currentLaunchVelocities, case2_currentLaunchAzimuth )                           \
                = singlePhasePlotDataExtraction( currentSolarPhase,
                                                 case2_solarPhase,
                                                 case2_launchVelocity,
                                                 case2_launchAzimuth )

        ( case3_currentLaunchVelocities, case3_currentLaunchAzimuth )                           \
                = singlePhasePlotDataExtraction( currentSolarPhase,
                                                 case3_solarPhase,
                                                 case3_launchVelocity,
                                                 case3_launchAzimuth )

        ( case4_currentLaunchVelocities, case4_currentLaunchAzimuth )                           \
                = singlePhasePlotDataExtraction( currentSolarPhase,
                                                 case4_solarPhase,
                                                 case4_launchVelocity,
                                                 case4_launchAzimuth )

        ( case5_currentLaunchVelocities, case5_currentLaunchAzimuth )                           \
                = singlePhasePlotDataExtraction( currentSolarPhase,
                                                 case5_solarPhase,
                                                 case5_launchVelocity,
                                                 case5_launchAzimuth )

        ( case6_currentLaunchVelocities, case6_currentLaunchAzimuth )                           \
                = singlePhasePlotDataExtraction( currentSolarPhase,
                                                 case6_solarPhase,
                                                 case6_launchVelocity,
                                                 case6_launchAzimuth )

        currentPlotHandle.plot( case1_currentLaunchAzimuth, case1_currentLaunchVelocities,
                            linestyle='--', marker='o', color='blue',
                            label=label1 if index == 0 else "__nolegend__" )
        currentPlotHandle.plot( case2_currentLaunchAzimuth, case2_currentLaunchVelocities,
                            linestyle='--', marker='o', color='red',
                            label=label2 if index == 0 else "__nolegend__" )
        currentPlotHandle.plot( case3_currentLaunchAzimuth, case3_currentLaunchVelocities,
                            linestyle='--', marker='o', color='green',
                            label=label3 if index == 0 else "__nolegend__" )
        currentPlotHandle.plot( case4_currentLaunchAzimuth, case4_currentLaunchVelocities,
                            linestyle='--', marker='o', color=colors.cnames['purple'],
                            label=label4 if index == 0 else "__nolegend__" )
        currentPlotHandle.plot( case5_currentLaunchAzimuth, case5_currentLaunchVelocities,
                            linestyle='--', marker='o', color='black',
                            label=label5 if index == 0 else "__nolegend__" )
        currentPlotHandle.plot( case6_currentLaunchAzimuth, case6_currentLaunchVelocities,
                            linestyle='--', marker='o', color=colors.cnames['orange'],
                            label=label6 if index == 0 else "__nolegend__" )

        currentPlotHandle.grid(True)
        currentPlotHandle.set_title( 'Solar phase = ' + str(currentSolarPhase) + ' [deg]' )
        currentPlotHandle.set_xlabel( 'Launch Azimuth [deg]' )
        currentPlotHandle.set_ylabel( 'Launch Velocity [m/s]' )
        # currentPlotHandle.xaxis.set_ticks( np.arange( 0.0, 360.0, 10.0 ) )
        currentPlotHandle.yaxis.set_ticks( np.arange( 0.0, 16.0, 1.0 ) )
        if index == 0:
            currentPlotHandle.legend( ).draggable( )

# Start timer.
start_time = time.time( )

## Plot the ellipsoidal shape of the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593

print "Extracting data from all databases now..."

# Connect to SQLite database for the perturbing case
try:
        database1 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                    + "multiple_launch_velocity_with_perturbations/"
                                    + "simulation_time_9_months/"
                                    + "3.2Density_1cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case1_solarPhase, case1_launchVelocity, case1_launchAzimuth ) = extractDataFromSQL( database1 )

ratio = areaToMassRatio( 1.0e-2, 3.2e3 )
label1 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

try:
        database2 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                    + "multiple_launch_velocity_with_perturbations/"
                                    + "simulation_time_9_months/"
                                    + "7.5Density_1cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case2_solarPhase, case2_launchVelocity, case2_launchAzimuth ) = extractDataFromSQL( database2 )

ratio = areaToMassRatio( 1.0e-2, 7.5e3 )
label2 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

try:
        database3 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                    + "multiple_launch_velocity_with_perturbations/"
                                    + "simulation_time_9_months/"
                                    + "3.2Density_5cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case3_solarPhase, case3_launchVelocity, case3_launchAzimuth ) = extractDataFromSQL( database3 )

ratio = areaToMassRatio( 5.0e-2, 3.2e3 )
label3 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

try:
        database4 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                    + "multiple_launch_velocity_with_perturbations/"
                                    + "simulation_time_9_months/"
                                    + "7.5Density_5cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case4_solarPhase, case4_launchVelocity, case4_launchAzimuth ) = extractDataFromSQL( database4 )

ratio = areaToMassRatio( 5.0e-2, 7.5e3 )
label4 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

try:
        database5 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                    + "multiple_launch_velocity_with_perturbations/"
                                    + "simulation_time_9_months/"
                                    + "3.2Density_10cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case5_solarPhase, case5_launchVelocity, case5_launchAzimuth ) = extractDataFromSQL( database5 )

ratio = areaToMassRatio( 10.0e-2, 3.2e3 )
label5 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

try:
        database6 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                    + "multiple_launch_velocity_with_perturbations/"
                                    + "simulation_time_9_months/"
                                    + "7.5Density_10cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case6_solarPhase, case6_launchVelocity, case6_launchAzimuth ) = extractDataFromSQL( database6 )

ratio = areaToMassRatio( 10.0e-2, 7.5e3 )
label6 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

print "Processing all data now...\n"

plotDataAndTable( case1_solarPhase, case1_launchVelocity, case1_launchAzimuth )
plotDataAndTable( case2_solarPhase, case2_launchVelocity, case2_launchAzimuth )
plotDataAndTable( case3_solarPhase, case3_launchVelocity, case3_launchAzimuth )
plotDataAndTable( case4_solarPhase, case4_launchVelocity, case4_launchAzimuth )
plotDataAndTable( case5_solarPhase, case5_launchVelocity, case5_launchAzimuth )
plotDataAndTable( case6_solarPhase, case6_launchVelocity, case6_launchAzimuth )

plotCompartiveCaptureCases( case1_solarPhase, case1_launchVelocity, case1_launchAzimuth, label1,
                            case2_solarPhase, case2_launchVelocity, case2_launchAzimuth, label2,
                            case3_solarPhase, case3_launchVelocity, case3_launchAzimuth, label3,
                            case4_solarPhase, case4_launchVelocity, case4_launchAzimuth, label4,
                            case5_solarPhase, case5_launchVelocity, case5_launchAzimuth, label5,
                            case6_solarPhase, case6_launchVelocity, case6_launchAzimuth, label6 )

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

