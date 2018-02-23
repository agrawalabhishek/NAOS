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

print "Extracting data now...\n"

data1 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude ),                        \
                                ROUND( initial_solar_phase_angle ),                         \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 );",                                      \
                     database )

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
                     database )

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
                     database )

data3.columns = [ 'vel_mag',                                                                \
                  'solar_phase_angle',                                                      \
                  'launch_azimuth' ]

capture_velocity_magnitude              = data3[ 'vel_mag' ]
capture_azimuth                         = data3[ 'launch_azimuth' ]
capture_solarPhase                      = data3[ 'solar_phase_angle' ]

if database:
    database.close( )

print "Capture cases = " + str( len( capture_velocity_magnitude ) )
print "Escape cases = " + str( len( escape_velocity_magnitude ) )
print "Reimpact cases = " + str( len( crash_velocity_magnitude ) )

print "Processing data now...\n"

## plot the histogram for number of particles escaped and reimpacted versus
## initial launch velocity
fig = plt.figure( )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )

plotHandles = [ ax1, ax2, ax3, ax4 ]

all_distinct_velocities = 16

uniqueSolarPhase = np.unique( escape_solarPhase )
for index in range( 0, len( uniqueSolarPhase ) ):
    currentPlotHandle = plotHandles[ index ]
    currentSolarPhase = uniqueSolarPhase[ index ]

    currentEscapeDataIndices = np.where( escape_solarPhase == currentSolarPhase )
    currentEscapeDataIndices = currentEscapeDataIndices[ 0 ]

    currentCrashDataIndices = np.where( crash_solarPhase == currentSolarPhase )
    currentCrashDataIndices = currentCrashDataIndices[ 0 ]

    currentCaptureDataIndices = np.where( capture_solarPhase == currentSolarPhase )
    currentCaptureDataIndices = currentCaptureDataIndices[ 0 ]

    escape_velocity_magnitude_plot = escape_velocity_magnitude[ currentEscapeDataIndices ]
    crash_velocity_magnitude_plot  = crash_velocity_magnitude[ currentCrashDataIndices ]
    capture_velocity_magnitude_plot = capture_velocity_magnitude[ currentCaptureDataIndices ]

    currentPlotHandle.hist( [crash_velocity_magnitude_plot, escape_velocity_magnitude_plot, capture_velocity_magnitude_plot],
                             bins=range(1, all_distinct_velocities+2),
                             histtype='bar',
                             color=['orange', 'crimson', 'green'],
                             label=['reimpact', 'escape', 'capture'],
                             align='left', alpha=0.7 )

    currentPlotHandle.set_xlim( 0, all_distinct_velocities+1 )
    currentPlotHandle.xaxis.set_ticks( np.arange(0, all_distinct_velocities+1, 1) )
    currentPlotHandle.set_xlabel( '$V_{initial}$ [m/s]' )
    currentPlotHandle.set_ylabel( 'Number of particles' )
    currentPlotHandle.set_title( 'Solar phase angle = ' + str(currentSolarPhase) + '[deg]' )
    currentPlotHandle.legend( ).draggable( )
    currentPlotHandle.grid( )

if database:
    database.close( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.suptitle( 'Regolith final fate histogram, Ellipsoid longest edge' )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
