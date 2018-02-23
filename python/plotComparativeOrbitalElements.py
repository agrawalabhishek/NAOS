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
import math
from scipy.interpolate import griddata
from decimal import Decimal

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
from matplotlib import animation
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

# asteroid's characteristic constants
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593
mu = 876514

# Set up the figure
fig = plt.figure( )
gs = gridspec.GridSpec( 3, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )
ax5 = plt.subplot( gs[ 4 ] )
ax6 = plt.subplot( gs[ 5 ] )

# Connect to SQLite database.
try:
       database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity/"
                                   + "simulation_time_9_months/"
                                   + "longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data for the non-perturbing case now...\n"

data = pd.read_sql( "SELECT     sma,                                                \
                                eccentricity,                                       \
                                inclination,                                        \
                                raan,                                               \
                                aop,                                                \
                                ta,                                                 \
                                ROUND( initial_velocity_magnitude ),                \
                                ROUND( launch_azimuth ),                            \
                                time                                                \
                     FROM       regolith_trajectory_results                         \
                     WHERE      ROUND( launch_azimuth ) = 185.0                     \
                     AND        ROUND( initial_velocity_magnitude ) = 5.0;",        \
                     database )

data.columns = [ 'sma',                                                 \
                 'eccentricity',                                        \
                 'inclination',                                         \
                 'raan',                                                \
                 'aop',                                                 \
                 'ta',                                                  \
                 'velocity_magnitude',                                  \
                 'launch_azimuth',                                      \
                 'time' ]

sma                 = data[ 'sma' ]
eccentricity        = data[ 'eccentricity' ]
inclination         = data[ 'inclination' ]
raan                = data[ 'raan' ]
aop                 = data[ 'aop' ]
ta                  = data[ 'ta' ]
velocityMagnitude   = data[ 'velocity_magnitude' ]
launchAzimuth       = data[ 'launch_azimuth' ]
t                   = data[ 'time' ]

endIndex = len( t )

if database:
    database.close( )

try:
    database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                               + "multiple_launch_velocity_with_perturbations/"
                               + "simulation_time_9_months/"
                               + "3.2Density_1cmSize/longestEdgePerturbations.db")
    # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
    #                            + "multiple_launch_velocity/"
    #                            + "simulation_time_9_months/"
    #                            + "longestEdge_v2.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data for the perturbing case now...\n"

lowerTimeDays = 0.0
lowerTime = lowerTimeDays * 24.0 * 60.0 * 60.0

upperTime = t[endIndex-1]
upperTimeDays = upperTime / ( 24.0 * 60.0 * 60.0 )

data = pd.read_sql( "SELECT     ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth ),                                        \
                                ROUND( initial_solar_phase_angle ),                             \
                                solar_phase_angle,                                              \
                                sma,                                                            \
                                eccentricity,                                                   \
                                inclination,                                                    \
                                raan,                                                           \
                                aop,                                                            \
                                ta,                                                             \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( launch_azimuth ) = 185.0                                 \
                     AND        ROUND( initial_velocity_magnitude ) = 5.0                       \
                     AND        ROUND( initial_solar_phase_angle ) = 135.0                      \
                     AND        time >= " + str(lowerTime)
                                + " AND time <= " + str(upperTime) + " ;",                      \
                     database )

data.columns = [ 'initial_velocity_magnitude',                                                  \
                 'launch_azimuth',                                                              \
                 'initial_solar_phase_angle',                                                   \
                 'solar_phase_angle',                                                           \
                 'sma',                                                                         \
                 'eccentricity',                                                                \
                 'inclination',                                                                 \
                 'raan',                                                                        \
                 'aop',                                                                         \
                 'ta',                                                                          \
                 'time' ]

# data = pd.read_sql( "SELECT     sma,                                                \
#                                 eccentricity,                                       \
#                                 inclination,                                        \
#                                 raan,                                               \
#                                 aop,                                                \
#                                 ta,                                                 \
#                                 ROUND( initial_solar_phase_angle ),                 \
#                                 solar_phase_angle,                                  \
#                                 ROUND( initial_velocity_magnitude ),                \
#                                 ROUND( launch_azimuth ),                            \
#                                 time                                                \
#                      FROM       regolith_trajectory_results                         \
#                      WHERE      ROUND( launch_azimuth ) = 185.0                     \
#                      AND        ROUND( initial_velocity_magnitude ) = 5.0;",        \
#                      database )

# data.columns = [ 'sma',                                                                         \
#                  'eccentricity',                                                                \
#                  'inclination',                                                                 \
#                  'raan',                                                                        \
#                  'aop',                                                                         \
#                  'ta',                                                                          \
#                  'initial_solar_phase_angle',                                                   \
#                  'solar_phase_angle',                                                           \
#                  'initial_velocity_magnitude',                                                  \
#                  'launch_azimuth',                                                              \
#                  'time' ]

if database:
    database.close( )

initialVelocityMagnitude                    = data[ 'initial_velocity_magnitude' ]
launchAzimuth                               = data[ 'launch_azimuth' ]
solarPhase                                  = data[ 'initial_solar_phase_angle' ]
changingSolarPhase                          = data[ 'solar_phase_angle' ]
t_withSolarPerturbations                    = data[ 'time' ]

sma_withSolarPerturbations                  = data[ 'sma' ]
eccentricity_withSolarPerturbations         = data[ 'eccentricity' ]
inclination_withSolarPerturbations          = data[ 'inclination' ]
raan_withSolarPerturbations                 = data[ 'raan' ]
aop_withSolarPerturbations                  = data[ 'aop' ]
ta_withSolarPerturbations                   = data[ 'ta' ]

print "Processing data now...\n"

print "Solar phase at start of the segment = " + str( changingSolarPhase[ 0 ] ) + " [deg]"

plt.suptitle( '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]', fontsize=12 )

## time in days
tDays = t / ( 24.0 * 60.0 * 60.0 )
tDays_withSolarPerturbations = t_withSolarPerturbations / ( 24.0 * 60.0 * 60.0 )

## plot the osculating orbital elements
# sma plot
ax1.plot( tDays_withSolarPerturbations, sma_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax1.plot( tDays, sma, color=colors.cnames['orange'], label='No solar perturbations' )
ax1.grid( True )
ax1.set_xlabel('time [days]')
ax1.set_ylabel('Semi-major axis [m]')
ax1.legend( ).draggable( )

# eccentricity
ax2.plot( tDays_withSolarPerturbations, eccentricity_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax2.plot( tDays, eccentricity, color=colors.cnames['orange'], label='No solar perturbations' )
ax2.grid( True )
ax2.set_xlabel('time [days]')
ax2.set_ylabel('Eccentricity')
ax2.legend( ).draggable( )

# inclination
ax3.plot( tDays_withSolarPerturbations, inclination_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax3.plot( tDays, inclination, color=colors.cnames['orange'], label='No solar perturbations' )
ax3.grid( True )
ax3.set_xlabel('time [days]')
ax3.set_ylabel('Inclination [deg]')
ax3.legend( ).draggable( )

# RAAN plot
ax4.plot( tDays_withSolarPerturbations, raan_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax4.plot( tDays, raan, color=colors.cnames['orange'], label='No solar perturbations' )
ax4.grid( True )
ax4.set_xlabel('time [days]')
ax4.set_ylabel('Right Ascension of Ascending Node [deg]')
ax4.legend( ).draggable( )

# AOP plot
ax5.plot( tDays_withSolarPerturbations, aop_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax5.plot( tDays, aop, color=colors.cnames['orange'], label='No solar perturbations' )
ax5.grid( True )
ax5.set_xlabel('time [days]')
ax5.set_ylabel('Argument of Periapsis [deg]')
ax5.legend( ).draggable( )

# TA plot
ax6.plot( tDays_withSolarPerturbations, ta_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax6.plot( tDays, ta, color=colors.cnames['orange'], label='No solar perturbations' )
ax6.grid( True )
ax6.set_xlabel('time [days]')
ax6.set_ylabel('True Anomaly [deg]')
ax6.legend( ).draggable( )

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
