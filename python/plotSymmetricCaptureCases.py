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

def plotEllipse( semiMajor, semiMinor, angleRange, plotHandle ):
        r = semiMajor * semiMinor / np.sqrt( ( semiMinor * np.cos( angleRange ) )**2
                                           + ( semiMajor * np.sin( angleRange ) )**2 )
        x = r * np.cos( angleRange )
        y = r * np.sin( angleRange )
        plotHandle.plot( x, y )
        plotHandle.grid( )

# Start timer.
start_time = time.time( )

## Plot the ellipsoidal shape of the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593

# Connect to SQLite database.
try:
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "captureTest2.db")
                                   # + "3.2Density_1cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

lowerTime = 0.0
upperTime = 270.0 * 24.0 * 60.0 * 60.0
# upperTime = 2000000.0
solarPhase = 45.0

# first capture case
data = pd.read_sql( "SELECT     position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                ROUND( initial_velocity_magnitude ),                            \
                                inertial_position_x,                                            \
                                inertial_position_y,                                            \
                                inertial_position_z,                                            \
                                ROUND( launch_azimuth ),                                        \
                                ROUND( launch_declination ),                                    \
                                eccentricity,                                                   \
                                total_energy,                                                   \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( initial_solar_phase_angle ) = " + str(solarPhase) + "    \
                     AND        time >= " + str(lowerTime)
                                + " AND time <= " + str(upperTime) + "                          \
                     AND        ROUND( launch_azimuth ) = 15.0                                  \
                     AND        ROUND( initial_velocity_magnitude ) = 8.0;",                    \
                     database )

data.columns = [ 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'velocity_magnitude',                                  \
                 'inertial_x',                                          \
                 'inertial_y',                                          \
                 'inertial_z',                                          \
                 'launch_azimuth',                                      \
                 'launch_declination',                                  \
                 'eccentricity',                                        \
                 'total_energy',                                        \
                 'time' ]

x                   = data[ 'x' ]
y                   = data[ 'y' ]
z                   = data[ 'z' ]
velocityMagnitude   = data[ 'velocity_magnitude' ]
inertial_x          = data[ 'inertial_x' ]
inertial_y          = data[ 'inertial_y' ]
inertial_z          = data[ 'inertial_z' ]
launchAzimuth       = data[ 'launch_azimuth' ]
launchDeclination   = data[ 'launch_declination' ]
eccentricity        = data[ 'eccentricity' ]
totalEnergy         = data[ 'total_energy' ]
t                   = data[ 'time' ]

# second capture case
data2 = pd.read_sql( "SELECT     position_x,                                                     \
                                 position_y,                                                     \
                                 position_z,                                                     \
                                 ROUND( initial_velocity_magnitude ),                            \
                                 inertial_position_x,                                            \
                                 inertial_position_y,                                            \
                                 inertial_position_z,                                            \
                                 ROUND( launch_azimuth ),                                        \
                                 ROUND( launch_declination ),                                    \
                                 eccentricity,                                                   \
                                 total_energy,                                                   \
                                 time                                                            \
                      FROM       regolith_trajectory_results                                     \
                      WHERE      ROUND( initial_solar_phase_angle ) = " + str(solarPhase) + "    \
                      AND        time >= " + str(lowerTime)
                                 + " AND time <= " + str(upperTime) + "                          \
                      AND        ROUND( launch_azimuth ) = 165.0                                 \
                      AND        ROUND( initial_velocity_magnitude ) = 8.0;",                    \
                      database )

data2.columns = [ 'x',                                                   \
                  'y',                                                   \
                  'z',                                                   \
                  'velocity_magnitude',                                  \
                  'inertial_x',                                          \
                  'inertial_y',                                          \
                  'inertial_z',                                          \
                  'launch_azimuth',                                      \
                  'launch_declination',                                  \
                  'eccentricity',                                        \
                  'total_energy',                                        \
                  'time' ]


x2                   = data2[ 'x' ]
y2                   = data2[ 'y' ]
z2                   = data2[ 'z' ]
velocityMagnitude2   = data2[ 'velocity_magnitude' ]
inertial_x2          = data2[ 'inertial_x' ]
inertial_y2          = data2[ 'inertial_y' ]
inertial_z2          = data2[ 'inertial_z' ]
launchAzimuth2       = data2[ 'launch_azimuth' ]
launchDeclination2   = data2[ 'launch_declination' ]
eccentricity2        = data2[ 'eccentricity' ]
totalEnergy2         = data2[ 'total_energy' ]
t2                   = data2[ 'time' ]

if database:
    database.close( )

print "Processing data now...\n"

# eccentricity and energy check for capture cases
fig = plt.figure( )
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
ax1.plot( t, eccentricity, label="Azimuth = " + str(launchAzimuth[ 0 ]) + " [deg]" )
ax1.plot( t2, eccentricity2, label="Azimuth = " + str(launchAzimuth2[ 0 ]) + " [deg]" )

ax2.plot( t, totalEnergy, label="Azimuth = " + str(launchAzimuth[ 0 ]) + " [deg]" )
ax2.plot( t2, totalEnergy2, label="Azimuth = " + str(launchAzimuth2[ 0 ]) + " [deg]" )

ax1.grid(True)
ax1.set_ylabel('Eccentricity')

ax2.grid(True)
ax2.set_ylabel('Energy')

## Inertial Frame
## Set up the figure
fig = plt.figure( figsize=( 7, 10 ) )
plt.suptitle( 'Ellipsoid longest edge, Inertial frame trajectory\n'
               + '$V_{launch}$ = ' + str( velocityMagnitude[ 0 ] ) + ' [m/s], '                 \
               + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '                   \
               + 'Solar phase = ' + str( solarPhase ) + ' [deg]\n'                              \
               + 'Time = ' + str( lowerTime / (24.0 * 60.0 * 60.0) ) + ' to '                   \
               + str( upperTime / (24.0 * 60.0 * 60.0) ) + ' [days]',                           \
               fontsize=10 )
gs = gridspec.GridSpec( 3, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

## Common format parameters
theta = np.linspace(0, 2 * np.pi, 360)
ellipsoidColor  = colors.cnames["slategray"]
trajectoryColor = colors.cnames["purple"]
textColor       = colors.cnames["black"]
startColor      = colors.cnames["darkgreen"]
endColor        = colors.cnames["darkred"]

## indicate ending point
endIndex = np.size( x )
endIndex2 = np.size( x2 )

## XY Projection
plotEllipse( alpha, beta, theta, ax1 )

ax1.plot( inertial_x, inertial_y, linewidth=0.5,
            label="Azimuth = " + str(launchAzimuth[ 0 ]) + " [deg]" )
ax1.plot( inertial_x2, inertial_y2, linewidth=0.5,
            label="Azimuth = " + str(launchAzimuth2[ 0 ]) + " [deg]" )

## indicate starting point
ax1.text( inertial_x[0], inertial_y[0], 'start', size=12, color=startColor )
ax1.text( inertial_x2[0], inertial_y2[0], 'start', size=12, color=startColor )

## indicate ending point
ax1.text( inertial_x[endIndex-1], inertial_y[endIndex-1], 'end', size=12, color=endColor )
ax1.text( inertial_x2[endIndex2-1], inertial_y2[endIndex2-1], 'end', size=12, color=endColor )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax1.grid( True )
ax1.legend( ).draggable( )
ax1.axis('equal')

## YZ Projection
plotEllipse( beta, gamma, theta, ax2 )

ax2.plot( inertial_y, inertial_z, linewidth=0.5,
            label="Azimuth = " + str(launchAzimuth[ 0 ]) + " [deg]" )
ax2.plot( inertial_y2, inertial_z2, linewidth=0.5,
            label="Azimuth = " + str(launchAzimuth2[ 0 ]) + " [deg]" )

## indicate starting point
ax2.text( inertial_y[0], inertial_z[0], 'start', size=12, color=startColor )
ax2.text( inertial_y2[0], inertial_z2[0], 'start', size=12, color=startColor )

## indicate ending point
ax2.text( inertial_y[endIndex-1], inertial_z[endIndex-1], 'end', size=12, color=endColor )
ax2.text( inertial_y2[endIndex2-1], inertial_z2[endIndex2-1], 'end', size=12, color=endColor )

## format axis and title
ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax2.grid( True )
ax2.legend( ).draggable( )
ax2.axis('equal')

## XZ Projection
plotEllipse( alpha, gamma, theta, ax3 )

ax3.plot( inertial_x, inertial_z, linewidth=0.5,
            label="Azimuth = " + str(launchAzimuth[ 0 ]) + " [deg]" )
ax3.plot( inertial_x2, inertial_z2, linewidth=0.5,
            label="Azimuth = " + str(launchAzimuth2[ 0 ]) + " [deg]" )

## indicate starting point
ax3.text( inertial_x[0], inertial_z[0], 'start', size=12, color=startColor )
ax3.text( inertial_x2[0], inertial_z2[0], 'start', size=12, color=startColor )

## indicate ending point
ax3.text( inertial_x[endIndex-1], inertial_z[endIndex-1], 'end', size=12, color=endColor )
ax3.text( inertial_x2[endIndex2-1], inertial_z2[endIndex2-1], 'end', size=12, color=endColor )

## format axis and title
ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax3.grid( True )
ax3.legend( ).draggable( )
ax3.axis('equal')

## seperate figure for the inertial frame 3D trajectory
## plot the 3d trajectory in inertial frame
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax2 = plt.subplot( gs[ 0 ], projection = '3d' )
ax2.plot( inertial_x, inertial_y, inertial_z, zdir = 'z', linewidth=0.5 )
ax2.plot( inertial_x2, inertial_y2, inertial_z2, zdir = 'z', linewidth=0.5 )

ax2.text( inertial_x[0], inertial_y[0], inertial_z[0], 'start', size=10, zorder=1, color='black' )
ax2.text( inertial_x[endIndex-1], inertial_y[endIndex-1], inertial_z[endIndex-1],
          'end', size=10, zorder=1, color='black' )

ax2.text( inertial_x2[0], inertial_y2[0], inertial_z2[0], 'start', size=10, zorder=1, color='black' )
ax2.text( inertial_x2[endIndex2-1], inertial_y2[endIndex2-1], inertial_z2[endIndex2-1],
          'end', size=10, zorder=1, color='black' )

ax2.set_aspect(1)
ax2.grid(True)
ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.set_zlabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
# ax2.set_title( 'Ellipsoid longest edge, Inertial frame trajectory\n'
#                + '$V_{launch}$ = ' + str( velocityMagnitude[ 0 ] ) + ' [m/s], '                 \
#                + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '                   \
#                + 'Solar phase = ' + str( solarPhase ) + ' [deg]\n'                              \
#                + 'Time = ' + str( lowerTime ) + ' to '                                          \
#                + str( upperTime / (24.0 * 60.0 * 60.0) ) + ' [days]',                           \
#                fontsize=10 )

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
