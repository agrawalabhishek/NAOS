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

# Connect to SQLite database.
try:
       database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity/"
                                   + "simulation_time_9_months/"
                                   + "longestEdge_v2.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data for the non-perturbing case now...\n"

data = pd.read_sql( "SELECT     kinetic_energy,                                     \
                                potential_energy,                                   \
                                total_energy,                                       \
                                true_potential_energy,                              \
                                true_total_energy,                                  \
                                ROUND( initial_velocity_magnitude ),                \
                                ROUND( launch_azimuth ),                            \
                                time                                                \
                     FROM       regolith_trajectory_results                         \
                     WHERE      ROUND( launch_azimuth ) = 165.0                     \
                     AND        ROUND( initial_velocity_magnitude ) = 8.0;",        \
                     database )

data.columns = [ 'kinetic_energy',                                      \
                 'potential_energy',                                    \
                 'total_energy',                                        \
                 'true_potential_energy',                               \
                 'true_total_energy',                                   \
                 'velocity_magnitude',                                  \
                 'launch_azimuth',                                      \
                 'time' ]

kinetic_energy          = data[ 'kinetic_energy' ]
potential_energy        = data[ 'potential_energy' ]
total_energy            = data[ 'total_energy' ]
true_potential_energy   = data[ 'true_potential_energy' ]
true_total_energy       = data[ 'true_total_energy' ]
velocityMagnitude       = data[ 'velocity_magnitude' ]
launchAzimuth           = data[ 'launch_azimuth' ]
t                       = data[ 'time' ]

endIndex = len( t )

if database:
    database.close( )

try:
    database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                               + "multiple_launch_velocity_with_perturbations/"
                               + "simulation_time_9_months/"
                               + "3.2Density_1cmSize/longestEdgePerturbations.db")

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
                                kinetic_energy,                                                 \
                                potential_energy,                                               \
                                total_energy,                                                   \
                                true_potential_energy,                                          \
                                true_total_energy,                                              \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( launch_azimuth ) = 165.0                                 \
                     AND        ROUND( initial_velocity_magnitude ) = 8.0                       \
                     AND        ROUND( initial_solar_phase_angle ) = 45.0                       \
                     AND        time >= " + str(lowerTime)
                                + " AND time <= " + str(upperTime) + " ;",                      \
                     database )

if database:
    database.close( )

data.columns = [ 'initial_velocity_magnitude',                                                  \
                 'launch_azimuth',                                                              \
                 'initial_solar_phase_angle',                                                   \
                 'solar_phase_angle',                                                           \
                 'kinetic_energy',                                                              \
                 'potential_energy',                                                            \
                 'total_energy',                                                                \
                 'true_potential_energy',                                                       \
                 'true_total_energy',                                                           \
                 'time' ]

initialVelocityMagnitude                    = data[ 'initial_velocity_magnitude' ]
launchAzimuth                               = data[ 'launch_azimuth' ]
solarPhase                                  = data[ 'initial_solar_phase_angle' ]
changingSolarPhase                          = data[ 'solar_phase_angle' ]
t_withSolarPerturbations                    = data[ 'time' ]

kinetic_energy_withSolarPerturbations          = data[ 'kinetic_energy' ]
potential_energy_withSolarPerturbations        = data[ 'potential_energy' ]
total_energy_withSolarPerturbations            = data[ 'total_energy' ]
true_potential_energy_withSolarPerturbations   = data[ 'true_potential_energy' ]
true_total_energy_withSolarPerturbations       = data[ 'true_total_energy' ]

print "Processing data now...\n"

print "Solar phase at start of the segment = " + str( changingSolarPhase[ 0 ] ) + " [deg]"

plt.suptitle( '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]', fontsize=12 )

## time in hrs
tDays = t / ( 60.0 * 60.0 )
tDays_withSolarPerturbations = t_withSolarPerturbations / ( 60.0 * 60.0 )

## plot the different energies
# kinetic energy
ax1.plot( tDays_withSolarPerturbations, kinetic_energy_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax1.plot( tDays, kinetic_energy, color=colors.cnames['orange'], linestyle='dotted', label='No solar perturbations' )
ax1.grid( True )
ax1.set_xlabel('Time [hrs]')
ax1.set_ylabel('Kinetic energy [$m^2/s^2$]')
ax1.legend( ).draggable( )

# potential energy
ax2.plot( tDays_withSolarPerturbations, potential_energy_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax2.plot( tDays, potential_energy, color=colors.cnames['orange'], linestyle='dotted', label='No solar perturbations' )
ax2.grid( True )
ax2.set_xlabel('Time [hrs]')
ax2.set_ylabel('Potential energy [$m^2/s^2$]')
ax2.legend( ).draggable( )

# total energy
ax4.plot( tDays_withSolarPerturbations, total_energy_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax4.plot( tDays, total_energy, color=colors.cnames['orange'], linestyle='dotted', label='No solar perturbations' )
ax4.grid( True )
ax4.set_xlabel('Time [hrs]')
ax4.set_ylabel('Total energy [$m^2/s^2$]')
ax4.legend( ).draggable( )

# true potential energy
ax3.plot( tDays_withSolarPerturbations, true_potential_energy_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax3.plot( tDays, true_potential_energy, color=colors.cnames['orange'], linestyle='dotted', label='No solar perturbations' )
ax3.grid( True )
ax3.set_xlabel('Time [hrs]')
ax3.set_ylabel('True potential energy [$m^2/s^2$]')
ax3.legend( ).draggable( )

# true total energy
ax5.plot( tDays_withSolarPerturbations, true_total_energy_withSolarPerturbations, color=colors.cnames['purple'], label='SRP + STBE' )
ax5.plot( tDays, true_total_energy, color=colors.cnames['orange'], linestyle='dotted', label='No solar perturbations' )
ax5.grid( True )
ax5.set_xlabel('Time [hrs]')
ax5.set_ylabel('True total energy [$m^2/s^2$]')
ax5.legend( ).draggable( )

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
