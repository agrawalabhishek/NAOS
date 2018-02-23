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
mu = 876514

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

data = pd.read_sql( "SELECT     trajectory_id,                                                  \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth ),                                        \
                                ROUND( initial_solar_phase_angle ),                             \
                                solar_phase_angle,                                              \
                                position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                srp_x,                                                          \
                                srp_y,                                                          \
                                srp_z,                                                          \
                                solarTide_x,                                                    \
                                solarTide_y,                                                    \
                                solarTide_z,                                                    \
                                gravAcc_x,                                                      \
                                gravAcc_y,                                                      \
                                gravAcc_z,                                                      \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( launch_azimuth ) = 185.0                                 \
                     AND        ROUND( initial_velocity_magnitude ) = 5.0                       \
                     AND        ROUND( initial_solar_phase_angle ) = 135.0;",                   \
                     database )

if database:
    database.close( )

data.columns = [ 'trajectory_id',                                                               \
                 'initial_velocity_magnitude',                                                  \
                 'launch_azimuth',                                                              \
                 'initial_solar_phase_angle',                                                   \
                 'solar_phase_angle',                                                           \
                 'position_x',                                                                  \
                 'position_y',                                                                  \
                 'position_z',                                                                  \
                 'srp_x',                                                                       \
                 'srp_y',                                                                       \
                 'srp_z',                                                                       \
                 'solarTide_x',                                                                 \
                 'solarTide_y',                                                                 \
                 'solarTide_z',                                                                 \
                 'gravAcc_x',                                                                   \
                 'gravAcc_y',                                                                   \
                 'gravAcc_z',                                                                   \
                 'time' ]

trajectory_id       = data[ 'trajectory_id' ]
velocityMagnitude   = data[ 'initial_velocity_magnitude' ]
launchAzimuth       = data[ 'launch_azimuth' ]
solarPhase          = data[ 'initial_solar_phase_angle' ]
changingSolarPhase  = data[ 'solar_phase_angle' ]
position_x          = data[ 'position_x' ]
position_y          = data[ 'position_y' ]
position_z          = data[ 'position_z' ]
srp_x               = data[ 'srp_x' ]
srp_y               = data[ 'srp_y' ]
srp_z               = data[ 'srp_z' ]
solarTide_x         = data[ 'solarTide_x' ]
solarTide_y         = data[ 'solarTide_y' ]
solarTide_z         = data[ 'solarTide_z' ]
gravAcc_x           = data[ 'gravAcc_x' ]
gravAcc_y           = data[ 'gravAcc_y' ]
gravAcc_z           = data[ 'gravAcc_z' ]
t                   = data[ 'time' ]

print "Processing data now...\n"

## Plot osculating time period variation with simulation time
fig = plt.figure( )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )

plt.suptitle('Gravity and perturbing accelerations, Ellipsoid longest edge \n'
             + 'Launch velocity = ' + str( velocityMagnitude[ 0 ] ) + ' [m/s], '
             + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
             + 'Initial Solar phase angle = ' + str( solarPhase[ 0 ] ) + ' [deg]' )

# convert time to days
t = t / ( 60.0 * 60.0 * 24.0 )

## plot srp magnitude with time
srpMag = np.sqrt( srp_x**2 + srp_y**2 + srp_z**2 )
ax1.plot( t, srpMag )
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
ax1.set_title('Change in SRP acceleration magnitude with time')
ax1.set_xlabel( 'Time [days]' )
ax1.set_ylabel( 'SRP acceleration $[m/s^2]$' )
ax1.grid(True)

## plot solar tidal effec with time
solarTidalEffect = np.sqrt( solarTide_x**2 + solarTide_y**2 + solarTide_z**2 )
ax2.plot( t, solarTidalEffect )
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
ax2.set_title('Change in Solar tidal acceleration magnitude with time')
ax2.set_xlabel( 'Time [days]' )
ax2.set_ylabel( 'Solar tidal acceleration $[m/s^2]$' )
ax2.grid(True)

## plot change in gravitational acceleration with time
gravAcc = np.sqrt( gravAcc_x**2 + gravAcc_y**2 + gravAcc_z**2 )
ax3.plot( t, gravAcc )
# ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax3.set_yscale( 'log' )
ax3.set_title('Change in Gravitational acceleration magnitude with time')
ax3.set_xlabel( 'Time [days]' )
ax3.set_ylabel( 'Gravitaitonal acceleration $[m/s^2]$' )
ax3.grid(True)

## calculate the net acceleration acting on the particle
net_x = srp_x + solarTide_x + gravAcc_x
net_y = srp_y + solarTide_y + gravAcc_y
net_z = srp_z + solarTide_z + gravAcc_z

netAcc = np.sqrt( net_x**2 + net_y**2 + net_z**2 )
ax4.plot( t, netAcc )
# ax4.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax4.set_yscale( 'log' )
ax4.set_title('Net acceleration acting on regolith')
ax4.set_xlabel( 'Time [days]' )
ax4.set_ylabel( 'Net acceleration $[m/s^2]$' )
ax4.grid(True)

## plot change in solar phase angle with time
fig = plt.figure(  )
gs = gridspec.GridSpec( 1, 1 )
ax5 = plt.subplot( gs[ 0 ] )
ax5.plot( t, changingSolarPhase )
# ax5.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax5.set_title('Evolution of Solar phase angle with time')
ax5.set_xlabel( 'Time [days]' )
ax5.set_ylabel( 'Solar phase angle [deg]' )
ax5.grid(True)

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
