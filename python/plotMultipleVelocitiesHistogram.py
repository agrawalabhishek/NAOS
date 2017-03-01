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
        # database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity_with_perturbations/phase_45/leadingEdge.db")
        database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

phaseAngle = 'N.A.'

## get data for escape cases
data1 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 );",                                      \
                     database )

data1.columns = [ 'vel_mag',                                                                \
                  'launch_azimuth' ]

escape_velocity_magnitude              = data1[ 'vel_mag' ]
escape_azimuth                         = data1[ 'launch_azimuth' ]

## get data for re-impact cases
data2 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( crash_flag = 1 );",                                       \
                     database )

data2.columns = [ 'vel_mag',                                                                \
                  'launch_azimuth' ]

crash_velocity_magnitude              = data2[ 'vel_mag' ]
crash_azimuth                         = data2[ 'launch_azimuth' ]

## get data for temporary capture cases
data3 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( end_flag = 1 );",                                         \
                     database )

data3.columns = [ 'vel_mag',                                                                \
                  'launch_azimuth' ]

capture_velocity_magnitude              = data3[ 'vel_mag' ]
capture_azimuth                         = data3[ 'launch_azimuth' ]

print "Capture cases = " + str( len( capture_velocity_magnitude ) )

## plot the histogram for number of particles escaped and reimpacted versus
## initial launch velocity
fig = plt.figure( )
ax3 = fig.add_subplot( 111 )

distinct_escape_velocities          = len( np.unique( escape_velocity_magnitude ) )
distinct_crash_velocities           = len( np.unique( crash_velocity_magnitude ) )
common_escape_crash_velocities      = len( np.intersect1d( escape_velocity_magnitude, crash_velocity_magnitude ) )

all_distinct_velocities = distinct_escape_velocities + distinct_crash_velocities - common_escape_crash_velocities

ax3.hist( [crash_velocity_magnitude, escape_velocity_magnitude],
          bins=range(1, all_distinct_velocities+2),
          histtype='bar',
          color=['orange', 'crimson'],
          label=['reimpact', 'escape'],
          align='left', alpha=0.7 )

ax3.set_xlim( 0, all_distinct_velocities+1 )
ax3.xaxis.set_ticks( np.arange(0, all_distinct_velocities+1, 1) )
ax3.set_xlabel( '$V_{initial}$ [m/s]' )
ax3.set_ylabel( 'Number of particles' )
ax3.set_title( 'Final fate histogram against launch velocities \n Spherical asteroid, Phase angle = ' + str(phaseAngle) + '[deg]' )
ax3.legend( ).draggable( )
ax3.grid( )

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
