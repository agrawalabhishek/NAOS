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

def plotEllipse( semiMajor, semiMinor, angleRange, plotHandle ):
        r = semiMajor * semiMinor / np.sqrt( ( semiMinor * np.cos( angleRange ) )**2
                                           + ( semiMajor * np.sin( angleRange ) )**2 )
        x = r * np.cos( angleRange )
        y = r * np.sin( angleRange )
        plotHandle.plot( x, y )
        plotHandle.grid( )

## Plot the ellipsoidal shape of the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593
mu = 876514

asteroidRotationPeriod = ( 2.0 * np.pi ) / Wz
asteroidRotationPeriod = asteroidRotationPeriod / ( 60.0 * 60.0 )

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

lowerTimeDays = 0.0
lowerTime = lowerTimeDays * 24.0 * 60.0 * 60.0

upperTimeDays = 5.0
upperTime = upperTimeDays * 24.0 * 60.0 * 60.0

data = pd.read_sql( "SELECT     trajectory_id,                                                  \
                                sma,                                                            \
                                eccentricity,                                                   \
                                inclination,                                                    \
                                raan,                                                           \
                                aop,                                                            \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth ),                                        \
                                ROUND( initial_solar_phase_angle ),                             \
                                solar_phase_angle,                                              \
                                position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                inertial_position_x,                                            \
                                inertial_position_y,                                            \
                                inertial_position_z,                                            \
                                srp_x,                                                          \
                                srp_y,                                                          \
                                srp_z,                                                          \
                                solarTide_x,                                                    \
                                solarTide_y,                                                    \
                                solarTide_z,                                                    \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( launch_azimuth ) = 185.0                                 \
                     AND        ROUND( initial_velocity_magnitude ) = 5.0                       \
                     AND        ROUND( initial_solar_phase_angle ) = 135.0                      \
                     AND        time >= " + str(lowerTime)
                                + " AND time <= " + str(upperTime) + " ;",                      \
                     database )

if database:
    database.close( )

data.columns = [ 'trajectory_id',                                                               \
                 'sma',                                                                         \
                 'eccentricity',                                                                \
                 'inclination',                                                                 \
                 'raan',                                                                        \
                 'aop',                                                                         \
                 'initial_velocity_magnitude',                                                  \
                 'launch_azimuth',                                                              \
                 'initial_solar_phase_angle',                                                   \
                 'solar_phase_angle',                                                           \
                 'position_x',                                                                  \
                 'position_y',                                                                  \
                 'position_z',                                                                  \
                 'inertial_position_x',                                                         \
                 'inertial_position_y',                                                         \
                 'inertial_position_z',                                                         \
                 'srp_x',                                                                       \
                 'srp_y',                                                                       \
                 'srp_z',                                                                       \
                 'solarTide_x',                                                                 \
                 'solarTide_y',                                                                 \
                 'solarTide_z',                                                                 \
                 'time' ]

trajectory_id       = data[ 'trajectory_id' ]
sma                 = data[ 'sma' ]
eccentricity        = data[ 'eccentricity' ]
inclination         = data[ 'inclination' ]
raan                = data[ 'raan' ]
aop                 = data[ 'aop' ]
velocityMagnitude   = data[ 'initial_velocity_magnitude' ]
launchAzimuth       = data[ 'launch_azimuth' ]
solarPhase          = data[ 'initial_solar_phase_angle' ]
changingSolarPhase  = data[ 'solar_phase_angle' ]
position_x          = data[ 'position_x' ]
position_y          = data[ 'position_y' ]
position_z          = data[ 'position_z' ]
inertial_position_x = data[ 'inertial_position_x' ]
inertial_position_y = data[ 'inertial_position_y' ]
inertial_position_z = data[ 'inertial_position_z' ]
srp_x               = data[ 'srp_x' ]
srp_y               = data[ 'srp_y' ]
srp_z               = data[ 'srp_z' ]
solarTide_x         = data[ 'solarTide_x' ]
solarTide_y         = data[ 'solarTide_y' ]
solarTide_z         = data[ 'solarTide_z' ]
t                   = data[ 'time' ]

print "Processing data now...\n"

print "Solar phase at start of the segment = " + str( changingSolarPhase[ 0 ] ) + " [deg]"

## Plot osculating time period variation with simulation time
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

timePeriod = ( 2.0 * np.pi ) * np.sqrt( sma**3 / mu )
timePeriod = timePeriod / ( 60.0 * 60.0 )
t = t / ( 60.0 * 60.0 * 24.0 )

ax1.plot( t, timePeriod,                                                                 \
          label='Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + '[deg], '
                + 'Launch Velocity = ' + str( velocityMagnitude[ 0 ] ) + '[m/s], '
                + 'Solar phase = ' + str( solarPhase[ 0 ] ) + '[deg]' )

ax1.axhline( y=asteroidRotationPeriod, color='red', label='Asteroid rotation period' )
ax1.axhline( y=5.0*asteroidRotationPeriod, color='orange', label='5:1 resonance' )
ax1.axhline( y=6.0*asteroidRotationPeriod, color='black', label='6:1 resonance' )
ax1.axhline( y=7.0*asteroidRotationPeriod, color=colors.cnames['purple'], label='7:1 resonance' )
ax1.axhline( y=8.0*asteroidRotationPeriod, color='green', label='8:1 resonance' )
ax1.axhline( y=9.0*asteroidRotationPeriod, color=colors.cnames['slategray'], label='9:1 resonance' )

ax1.set_ylabel('Osculating orbital time period [hrs]')
ax1.set_xlabel('Simulation Time [Days]')
ax1.set_yscale('log')
# ax1.set_xscale('log')
ax1.grid( True )
ax1.legend( ).draggable( )

## plot transition zone trajectory segment with srp and tidal vectors in body frame
fig = plt.figure( )
gs = gridspec.GridSpec( 3, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 2 ] )
ax3 = plt.subplot( gs[ 4 ] )
ax4 = plt.subplot( gs[ 1 ] )
ax5 = plt.subplot( gs[ 3 ] )
ax6 = plt.subplot( gs[ 5 ] )

plt.suptitle( 'Trajectory segment with direction of perturbing accelerations' )

theta = np.linspace(0, 2 * np.pi, 300)

plotEllipse( alpha, beta, theta, ax1 )
plotEllipse( beta, gamma, theta, ax2 )
plotEllipse( alpha, gamma, theta, ax3 )
plotEllipse( alpha, beta, theta, ax4 )
plotEllipse( beta, gamma, theta, ax5 )
plotEllipse( alpha, gamma, theta, ax6 )

ax1.plot( position_x, position_y, linewidth=1 )
ax1.scatter( position_x[ 0 ], position_y[ 0 ], marker='*', s=12, c='k' )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')

ax2.plot( position_y, position_z, linewidth=1 )
ax2.scatter( position_y[ 0 ], position_z[ 0 ], marker='*', s=12, c='k' )
ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')

ax3.plot( position_x, position_z, linewidth=1 )
ax3.scatter( position_x[ 0 ], position_z[ 0 ], marker='*', s=12, c='k' )
ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')

ax4.plot( position_x, position_y, linewidth=1 )
ax4.scatter( position_x[ 0 ], position_y[ 0 ], marker='*', s=12, c='k' )
ax4.set_xlabel('x [m]')
ax4.set_ylabel('y [m]')

ax5.plot( position_y, position_z, linewidth=1 )
ax5.scatter( position_y[ 0 ], position_z[ 0 ], marker='*', s=12, c='k' )
ax5.set_xlabel('y [m]')
ax5.set_ylabel('z [m]')

ax6.plot( position_x, position_z, linewidth=1 )
ax6.scatter( position_x[ 0 ], position_z[ 0 ], marker='*', s=12, c='k' )
ax6.set_xlabel('x [m]')
ax6.set_ylabel('z [m]')

# indices1 = range( 0, 20000, 500 )
# indices2 = range( 20000, len( position_x ), 10 )
# indices = indices1 + indices2
indices = range( 0, len( position_x ), 20 )
vector_loc_x = position_x[ indices ]
vector_loc_y = position_y[ indices ]
vector_loc_z = position_z[ indices ]
plot_srp_x = srp_x[ indices ]
plot_srp_y = srp_y[ indices ]
plot_srp_z = srp_z[ indices ]
plot_solarTide_x = solarTide_x[ indices ]
plot_solarTide_y = solarTide_y[ indices ]
plot_solarTide_z = solarTide_z[ indices ]

ax1.quiver( vector_loc_x, vector_loc_y,
            plot_srp_x, plot_srp_y,
            pivot='tail',
            color=colors.cnames[ 'red' ], linestyles='solid', lw=1, label='SRP' )


ax2.quiver( vector_loc_y, vector_loc_z,
            plot_srp_y, plot_srp_z,
            pivot='tail',
            color=colors.cnames[ 'red' ], linestyles='solid', lw=1, label='SRP' )


ax3.quiver( vector_loc_x, vector_loc_z,
            plot_srp_x, plot_srp_z,
            pivot='tail',
            color=colors.cnames[ 'red' ], linestyles='solid', lw=1, label='SRP' )

ax4.quiver( vector_loc_x, vector_loc_y,
            plot_solarTide_x, plot_solarTide_y,
            pivot='tail',
            color=colors.cnames[ 'orange' ], linestyles='solid', lw=1, label='Sun Third-Body effect' )

ax5.quiver( vector_loc_y, vector_loc_z,
            plot_solarTide_y, plot_solarTide_z,
            pivot='tail',
            color=colors.cnames[ 'orange' ], linestyles='solid', lw=1, label='Sun Third-Body effect' )

ax6.quiver( vector_loc_x, vector_loc_z,
            plot_solarTide_x, plot_solarTide_z,
            pivot='tail',
            color=colors.cnames[ 'orange' ], linestyles='solid', lw=1, label='Sun Third-Body effect' )

## plot transition zone trajectory segment in inertial frame
fig = plt.figure( )
gs = gridspec.GridSpec( 3, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

plt.suptitle( 'Transition zone trajectory segment (Inertial Frame)' )

plotEllipse( alpha, beta, theta, ax1 )
plotEllipse( beta, gamma, theta, ax2 )
plotEllipse( alpha, gamma, theta, ax3 )

ax1.plot( inertial_position_x, inertial_position_y, linewidth=0.5 )
ax1.scatter( inertial_position_x[ 0 ], inertial_position_y[ 0 ], marker='*', s=8, c='k' )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')

ax2.plot( inertial_position_y, inertial_position_z, linewidth=0.5 )
ax2.scatter( inertial_position_y[ 0 ], inertial_position_z[ 0 ], marker='*', s=8, c='k' )
ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')

ax3.plot( inertial_position_x, inertial_position_z, linewidth=0.5 )
ax3.scatter( inertial_position_x[ 0 ], inertial_position_z[ 0 ], marker='*', s=8, c='k' )
ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')

## plot osculating orbital elements now and range for the segment
fig = plt.figure( )
gs = gridspec.GridSpec( 3, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )
ax5 = plt.subplot( gs[ 4 ] )
ax6 = plt.subplot( gs[ 5 ] )

## sma plot
ax1.plot( t, sma )
ax1.set_ylabel('Semi-major axis [m]')
ax1.set_xlabel('Simulation Time [Days]')
ax1.set_yscale('log')
# ax1.set_xscale('log')
ax1.grid( True )
ax1.legend( ).draggable( )

## eccentricity plot
ax2.plot( t, eccentricity )
ax2.set_ylabel('Eccentricity')
ax2.set_xlabel('Simulation Time [Days]')
# ax2.set_yscale('log')
# ax2.set_xscale('log')
ax2.grid( True )
ax2.legend( ).draggable( )

## inclination plot
ax3.plot( t, inclination )
ax3.set_ylabel('Inclination [deg]')
ax3.set_xlabel('Simulation Time [Days]')
# ax3.set_yscale('log')
# ax3.set_xscale('log')
ax3.grid( True )
ax3.legend( ).draggable( )

## RAAN plot
ax4.plot( t, raan )
ax4.set_ylabel('RAAN [deg]')
ax4.set_xlabel('Simulation Time [Days]')
# ax4.set_yscale('log')
# ax4.set_xscale('log')
ax4.grid( True )
ax4.legend( ).draggable( )

## AOP plot
ax5.plot( t, aop )
ax5.set_ylabel('AOP [deg]')
ax5.set_xlabel('Simulation Time [Days]')
# ax5.set_yscale('log')
# ax5.set_xscale('log')
ax5.grid( True )
ax5.legend( ).draggable( )

## range plot for the trajectory segment
segmentRange = np.sqrt( position_x**2 + position_y**2 + position_z**2 )
ax6.plot( t, segmentRange )
ax6.set_ylabel( 'Altitude [m]' )
ax6.set_xlabel( 'Simulation Time [days]' )
ax6.set_yscale( 'log' )
ax6.grid( True )

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
