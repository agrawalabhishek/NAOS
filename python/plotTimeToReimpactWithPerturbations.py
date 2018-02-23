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

## Operations
# Connect to SQLite database.
try:
    database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                               + "multiple_launch_velocity_with_perturbations/"
                               + "simulation_time_9_months/"
                               + "longestEdge_3P2Density_1cmRadius.db")
                               # + "3.2Density_1cmSize/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

data = pd.read_sql( "SELECT    ROUND( initial_solar_phase_angle )                          \
                    FROM       regolith_trajectory_results                                 \
                    WHERE      ( crash_flag = 1 );",                                       \
                    database )

data.columns = [ 'initial_solar_phase_angle' ]

solarPhases = data[ 'initial_solar_phase_angle' ]
uniqueSolarPhases = np.unique( solarPhases )

## plot escape time
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

    data1 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude ),                        \
                                    time,                                                       \
                                    ROUND( launch_azimuth )                                     \
                         FROM       regolith_trajectory_results                                 \
                         WHERE      crash_flag = 1                                              \
                         AND        ROUND( initial_solar_phase_angle ) = "
                                    + str( currentSolarPhase ) + ";",                           \
                         database )

    data1.columns = [ 'vel_mag',                                                                \
                      'time',                                                                   \
                      'launch_azimuth' ]

    crash_velocity_magnitude              = data1[ 'vel_mag' ]
    crash_time                            = data1[ 'time' ]
    crash_azimuth                         = data1[ 'launch_azimuth' ]

    print "Processing data...\n"

    currentUniqueVelocities = np.unique( crash_velocity_magnitude )
    for newIndex in range( 0, len( currentUniqueVelocities ) ):
        currentVelocity = currentUniqueVelocities[ newIndex ]

        newCrashDataIndices = np.where( crash_velocity_magnitude == currentVelocity )
        newCrashDataIndices = newCrashDataIndices[ 0 ]

        plotAzimuths = crash_azimuth[ newCrashDataIndices ]
        plotTime = crash_time[ newCrashDataIndices ]

        currentPlotHandle.scatter( plotAzimuths, plotTime, s=5, c=colors[newIndex],             \
                                   edgecolors='face',                                           \
                                   label='$V_{Launch}$ = ' + str( currentVelocity ) + ' [m/s]' )

        currentPlotHandle.xaxis.set_ticks( np.arange(0, 361, 30) )
        currentPlotHandle.set_xlabel( 'Launch azimuth' )
        currentPlotHandle.set_ylabel( 'Time [s]' )
        currentPlotHandle.set_yscale( 'log' )
        currentPlotHandle.set_title( 'Solar phase angle = ' + str(currentSolarPhase) + '[deg]' )
        currentPlotHandle.legend( markerscale=7 ).draggable( )
        currentPlotHandle.grid( True )

if database:
    database.close( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.suptitle( 'Time for regolith to reimpact, Ellipsoid Longest edge' )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
