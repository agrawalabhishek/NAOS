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

## Operations
# Connect to SQLite database.
try:
    database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                               + "multiple_launch_velocity_with_perturbations/"
                               + "simulation_time_9_months/"
                               + "3.2Density_1cmSize/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

data = pd.read_sql( "SELECT    ROUND( initial_solar_phase_angle ),                          \
                               trajectory_id                                                \
                    FROM       regolith_trajectory_results                                  \
                    WHERE      ( crash_flag = 1 );",                                        \
                    database )

data.columns = [ 'initial_solar_phase_angle',
                 'trajectory_id' ]

solarPhases = data[ 'initial_solar_phase_angle' ]
trajectory_id = data[ 'trajectory_id' ]
trajectory_id = tuple( trajectory_id.tolist( ) )

uniqueSolarPhases = np.unique( solarPhases )

fig = plt.figure( )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )

plotHandles = [ ax1, ax2, ax3, ax4 ]

colors = plt.cm.Vega20( np.linspace( 0, 1, 16 ) )

for solarIndex in range( 0, len( uniqueSolarPhases ) ):
    currentPlotHandle = plotHandles[ solarIndex ]
    currentSolarPhase = uniqueSolarPhases[ solarIndex ]

    data1 = pd.read_sql( "SELECT    trajectory_id,                                              \
                                    position_x,                                                 \
                                    position_y,                                                 \
                                    position_z,                                                 \
                                    ROUND( initial_velocity_magnitude ),                        \
                                    ROUND( launch_azimuth )                                     \
                         FROM       regolith_trajectory_results                                 \
                         WHERE      trajectory_id IN " + str( trajectory_id ) + "               \
                         AND        ROUND( initial_solar_phase_angle ) = "
                                    + str( currentSolarPhase ) + ";",                           \
                         database )

    data1.columns = [ 'trajectory_id',                                       \
                      'x',                                                   \
                      'y',                                                   \
                      'z',                                                   \
                      'init_vel_mag',                                        \
                      'launch_azimuth' ]

    print "Processing data for solar phase angle = " + str( currentSolarPhase ) + " [deg]\n"

    x                               = data1[ 'x' ]
    y                               = data1[ 'y' ]
    z                               = data1[ 'z' ]
    initial_velocity                = data1[ 'init_vel_mag' ]
    azimuth                         = data1[ 'launch_azimuth' ]
    currentSolarPhaseTrajectoryIds  = data1[ 'trajectory_id' ]

    radialDistances = np.sqrt( x**2 + y**2 + z**2 )

    unique_reimpact_velocities = np.unique( initial_velocity )
    for index in range( 0, len( unique_reimpact_velocities ) ):
        current_velocity = unique_reimpact_velocities[ index ]
        current_velocity_indices = np.where( initial_velocity == current_velocity )
        current_velocity_indices = current_velocity_indices[ 0 ]

        current_velocity_azimuths = azimuth[ current_velocity_indices ]
        current_velocity_radial_distances = radialDistances[ current_velocity_indices ]

        current_velocity_unique_azimuths = np.unique( current_velocity_azimuths )
        for index2 in range( 0, len( current_velocity_unique_azimuths ) ):
            current_azimuth = current_velocity_unique_azimuths[ index2 ]
            dataIndices = np.where( (azimuth == current_azimuth) & (initial_velocity == current_velocity) )
            dataIndices = dataIndices[ 0 ]
            currentPlotHandle.scatter( current_azimuth, max(radialDistances[dataIndices]),
                                       s=5, c=colors[ index ],
                                       edgecolors='face',
                                       label='$V_{launch}$ = ' + str( current_velocity ) + ' [m/s]' if index2 == 0 else '_nolegend_' )

            currentPlotHandle.legend( markerscale=2 ).draggable( )

    currentPlotHandle.set_xlabel('Launch azimuth [degree]')
    currentPlotHandle.set_ylabel('Max. radial distance [m]')
    currentPlotHandle.set_yscale('log')
    # cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )
    # currentPlotHandle.set_xlim( -180.0, 180.0 )
    # currentPlotHandle.set_yticks( np.arange( -90.0, 90.0, 10.0 ) )
    currentPlotHandle.set_title( 'Solar phase angle = ' + str( currentSolarPhase ) + ' [deg]' )
    currentPlotHandle.grid( True )

if database:
    database.close( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.suptitle( 'Max. radial distance in orbit for re-impact trajectories \n Ellipsoid longest edge' )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
