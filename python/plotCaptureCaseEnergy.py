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
    # get data for temporary capture cases
    data3 = pd.read_sql( "SELECT    trajectory_id,                                              \
                                    ROUND( initial_velocity_magnitude ),                        \
                                    ROUND( initial_solar_phase_angle ),                         \
                                    ROUND( launch_azimuth )                                     \
                         FROM       regolith_trajectory_results                                 \
                         WHERE      ( end_flag = 1 );",                                         \
                         databaseConnect )

    data3.columns = [ 'trajectory_id',                                                          \
                      'vel_mag',                                                                \
                      'solar_phase_angle',                                                      \
                      'launch_azimuth' ]

    uniquecaptureTrajectoryIds              = data3[ 'trajectory_id' ]
    capture_velocity_magnitude              = data3[ 'vel_mag' ]
    capture_azimuth                         = data3[ 'launch_azimuth' ]
    capture_solarPhase                      = data3[ 'solar_phase_angle' ]

    uniqueCaptureTrajectoryIds = tuple( uniquecaptureTrajectoryIds.tolist( ) )

    data4 = pd.read_sql( "SELECT    trajectory_id,                                              \
                                    total_energy,                                               \
                                    time                                                        \
                          FROM      regolith_trajectory_results                                 \
                          WHERE     trajectory_id IN " + str( uniqueCaptureTrajectoryIds ) + ";",\
                          databaseConnect )

    data4.columns = [ 'trajectory_id',                                                          \
                      'total_energy',                                                           \
                      'time' ]

    capture_total_energy        = data4[ 'total_energy' ]
    capture_time                = data4[ 'time' ]
    capture_trajectoryID        = data4[ 'trajectory_id' ]

    if databaseConnect:
        databaseConnect.close( )

    return ( capture_velocity_magnitude, capture_azimuth, capture_solarPhase, uniqueCaptureTrajectoryIds,
             capture_total_energy, capture_time, capture_trajectoryID )

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

## Operations
# Connect to SQLite database.
print "Extracting data now...\n"

try:
    database1 = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                + "multiple_launch_velocity_with_perturbations/"
                                + "simulation_time_9_months/"
                                + "3.2Density_1cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

( case1_launchVel, case1_azimuth, case1_solarPhase, case1_uniqueTrajectoryIds,
    case1_total_energy, case1_time, case1_trajectoryId ) = extractDataFromSQL( database1 )
ratio = areaToMassRatio( 1.0e-2, 3.2e3 )
label1 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'

####################################################################################################
# plot the total energy for number of particles CAPTURED against total flight time
fig = plt.figure( figsize=(12, 6) )
gs = gridspec.GridSpec( 1, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )

# for metadata
# ax3 = plt.subplot( gs[ 2 ] )
# ax3.axis( 'off' )
# metadata_table = []
# columnLabels = [ "Index", "$V_{launch}$ [m/s]", "Azimuth [deg]", "Solar phase [deg]" ]

for index in range( 0, len( case1_uniqueTrajectoryIds ) ):
    currentTrajectory = case1_uniqueTrajectoryIds[ index ]
    dataIndices = np.where( case1_trajectoryId == currentTrajectory )
    dataIndices = dataIndices[ 0 ]

    currentLaunchVelocity = case1_launchVel[ index ]
    currentLaunchAzimuth = case1_azimuth[ index ]
    currentSolarPhase = case1_solarPhase[ index ]

    plotEnergyValues = case1_total_energy[ dataIndices ]
    plotEnergyValues = plotEnergyValues.tolist( )

    plotTimeValues = case1_time[ dataIndices ] / ( 24.0 * 60.0 * 60.0 )
    plotTimeValues = plotTimeValues.tolist( )

    # relativeEnergy = []
    # plotTime = []
    # for energyIndex in range( 0, len( plotEnergyValues )-1 ):
    #     diff = plotEnergyValues[ energyIndex + 1 ] - plotEnergyValues[ energyIndex ]
    #     relativeEnergy.append( diff )
    #     plotTime.append( plotTimeValues[ energyIndex + 1 ] )

    # ax1.plot( plotTimeValues, plotEnergyValues,
    #           label='$V_{launch}$ [m/s] = ' + str(currentLaunchVelocity) +
    #                 'Azimuth [deg] = ' + str(currentLaunchAzimuth) +
    #                 'Solar phase [deg] = ' + str(currentSolarPhase) )
    plotHandle = ax1.plot( plotTimeValues, plotEnergyValues, label=str(index+1) )

    # metadata_table.append( [ index+1,                                                           \
    #                          currentLaunchVelocity,                                             \
    #                          currentLaunchAzimuth,                                              \
    #                          currentSolarPhase ] )

    ax1.set_xlabel( 'Time [days]' )
    ax1.set_ylabel( 'Total Energy [$m^2/s^2$]' )
    ax1.set_title( label1 )
    ax1.legend( ).draggable( )
    ax1.grid(True)

    ax2.plot( plotTimeValues, plotEnergyValues )
    ax2.set_xlabel( 'Time [days]' )
    ax2.set_ylabel( 'Total Energy [$m^2/s^2$]' )
    ax2.set_title( label1 )
    ax2.grid(True)

# table = ax3.table( cellText = metadata_table,
#                    colLabels = columnLabels,
#                    cellLoc = 'center',
#                    loc = 'center' )
# table.auto_set_font_size(False)
# table.set_fontsize( 10 )
# table_properties = table.properties( )
# table_cells = table_properties[ 'child_artists' ]
# for cell in table_cells: cell.set_height( 0.10 )
# cell_dict = table.get_celld( )

plt.suptitle( 'Capture case energy plots, Ellipsoid longest edge' )

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
