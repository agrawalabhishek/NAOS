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
        # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
        #                            + "multiple_launch_velocity_with_perturbations/"
        #                            + "simulation_time_9_months/"
        #                            + "3.2Density_1cmRadius/longestEdgePerturbations.db")
        # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
        #                            + "multiple_launch_velocity_with_perturbations/"
        #                            + "simulation_time_9_months/"
        #                            + "test.db")
        database = sqlite3.connect("../data/guarantee_escape_speed/"
                                   + "longest_edge/"
                                   + "normalLaunch.db")
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/"
        #                             + "spherical_asteroid/longestEdge.db" )

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

lowerTime = 0.0
upperTime = 270.0 * 24.0 * 60.0 * 60.0
solarPhase = 315.0

data = pd.read_sql( "SELECT     position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                ROUND( initial_velocity_magnitude ),                            \
                                inertial_position_x,                                            \
                                inertial_position_y,                                            \
                                inertial_position_z,                                            \
                                ROUND( launch_azimuth ),                                        \
                                ROUND( launch_declination ),                                    \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( initial_solar_phase_angle ) = " + str(solarPhase) + "    \
                     AND        time >= " + str(lowerTime)
                                + " AND time <= " + str(upperTime) + "                          \
                     AND        ROUND( launch_declination ) = 0.0                               \
                     AND        ROUND( initial_velocity_magnitude ) = 9.0;",                    \
                     database )

# data = pd.read_sql( "SELECT     position_x,                                                     \
#                                 position_y,                                                     \
#                                 position_z,                                                     \
#                                 ROUND( initial_velocity_magnitude ),                            \
#                                 inertial_position_x,                                            \
#                                 inertial_position_y,                                            \
#                                 inertial_position_z,                                            \
#                                 ROUND( launch_azimuth ),                                        \
#                                 ROUND( launch_declination ),                                    \
#                                 time                                                            \
#                      FROM       regolith_trajectory_results                                     \
#                      WHERE      ROUND( launch_azimuth ) = 246.0                                 \
#                      AND        ROUND( initial_velocity_magnitude ) = 11.0;",                   \
#                      database )

data.columns = [ 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'velocity_magnitude',                                  \
                 'inertial_x',                                          \
                 'inertial_y',                                          \
                 'inertial_z',                                          \
                 'launch_azimuth',                                      \
                 'launch_declination',                                  \
                 'time' ]

if database:
    database.close( )

x                   = data[ 'x' ]
y                   = data[ 'y' ]
z                   = data[ 'z' ]
velocityMagnitude   = data[ 'velocity_magnitude' ]
inertial_x          = data[ 'inertial_x' ]
inertial_y          = data[ 'inertial_y' ]
inertial_z          = data[ 'inertial_z' ]
launchAzimuth       = data[ 'launch_azimuth' ]
launchDeclination   = data[ 'launch_declination' ]
t                   = data[ 'time' ]

print "Processing data now...\n"

## sanity check for the initial position coordinate
ellipsoidSolution = ( np.square(x[0]) / alpha**2 ) + ( np.square(y[0]) / beta**2 ) + ( np.square(z[0]) / gamma**2 )
ellipsoidSolution = ellipsoidSolution - 1.0
if ellipsoidSolution > 1.0e-12:
    print "launch point not on ellipsoid"

## sanity check to ensure particle never impacted the surface and continued on with its trajectory
xSquare = x**2
ySquare = y**2
zSquare = z**2
impactCheck = (xSquare / alpha**2) + (ySquare / beta**2) + (zSquare / gamma**2) - 1.0

rangeCheck = np.sqrt( x**2 + y**2 + z**2 )

# firstRun = True
# impactCheckFlag = False
# for index in range( 1, len( impactCheck ) ):
#     if impactCheck[index] <= 0.0:
#         if firstRun == True:
#             fig = plt.figure( )
#             ax1 = plt.subplot(211)
#             ax2 = plt.subplot(212)
#             theta = np.linspace(0, 2 * np.pi, 360)
#             plotEllipse( alpha, beta, theta, ax1 )
#             ax1.plot( x, y )
#             ax2.plot( t, rangeCheck )
#             firstRun = False
#             print "\nParticle impacted yet trajectory continued in simulation - This is an error!!\n"
#             print "Impact at time = " + str( t[index] ) + " [s]\n"

#         ax1.scatter( x[index], y[index], c='red' )
#         ax2.scatter( t[index], rangeCheck[index], c='red' )
#         impactCheckFlag = True

# if impactCheckFlag == True:
#     plt.show( )
#     sys.exit( )

fig = plt.figure( )
ax1 = plt.subplot(111)
ax1.plot( t, rangeCheck )
ax1.axhline( beta, color='red' )
ax1.grid( True )
ax1.set_xlabel( "Time [s]" )
ax1.set_ylabel( "Range [m]" )
ax1.set_title( "Range to particle v/s time in orbit" )

## Set up the figure for the main plots
fig = plt.figure( figsize=( 6, 6 ) )
# ax1 = fig.add_subplot( 211, projection = '3d' )
# ax2 = fig.add_subplot( 212, frameon = False )
gs = gridspec.GridSpec( 2, 1, height_ratios = [ 1, 1 ] )
ax1 = plt.subplot( gs[ 0 ], projection = '3d' )
ax2 = plt.subplot( gs[ 1 ], projection = '3d' )

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

newColor = colors.cnames["brown"]
surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         rstride=5, cstride=5, alpha=1.0,
                         color=newColor )
# surf.set_facecolor( newColor )
surf.set_linewidth( 1.0 )

surf2 = ax2.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         rstride=5, cstride=5, alpha=1.0,
                         color=newColor )
# surf2.set_facecolor( newColor )
surf2.set_linewidth( 1.0 )

## Plot 3D trajectory of the orbiting particle
ax1.plot( x, y, z, zdir = 'z', color=colors.cnames["purple"], linewidth=1.0 )

# velocity vector
# ax1.quiver3D( x[10000], y[10000], z[10000],
#               vx[10000], vy[10000], vz[10000],
#               length=10000.0, lw = 1, pivot='tail', arrow_length_ratio=0.2,
#               color=colors.cnames["red"], linestyles='solid' )

# a position vector
# radius = np.sqrt( x[10000]**2 + y[10000]**2 + z[10000]**2 )
# ax1.quiver3D( 0.0, 0.0, 0.0,
#               x[10000], y[10000], z[10000],
#               length=radius, lw = 1, pivot='tail', arrow_length_ratio=0.2,
#               color=colors.cnames['green'], linestyles='solid' )

## indicate starting point
# ax1.scatter( x[0], y[0], z[0], 'g^' )
ax1.text( x[0], y[0], z[0], 'start', size=10, zorder=1, color=colors.cnames["black"] )

## indicate ending point
endIndex = np.size( x )
# ax1.scatter( x[endIndex - 1], y[endIndex - 1], z[endIndex - 1], 'k^' )
ax1.text( x[endIndex-1], y[endIndex-1], z[endIndex-1],
          'end', size=10, zorder=1,
          color=colors.cnames["black"] )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)

# ax1.set_title( 'Particle trajectory around asteroid Eros (Body frame) \n'
#                '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '                 \
#                'Launch azimuth=' + str( launchAzimuth[ 0 ] ) + '[deg], '                    \
#                'time=' + str( t[ endIndex-1 ] / (24.0*60.0*60.0) ) + '[day(s)]' )

# ax1.set_title( 'Ellipsoid longest edge, Rotating frame trajectory\n'
#                + '$V_{launch}$ = ' + str( velocityMagnitude[ 0 ] ) + ' [m/s], '                 \
#                + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '                   \
#                + 'Solar phase = ' + str( solarPhase ) + ' [deg]\n'                              \
#                + 'Time = ' + str( lowerTime ) + ' to '                                          \
#                + str( upperTime / (24.0 * 60.0 * 60.0) ) + ' [days]',                           \
#                fontsize=10 )

ax1.set_title( 'Ellipsoid longest edge, Rotating frame trajectory\n'
               + '$V_{launch}$ = ' + str( velocityMagnitude[ 0 ] ) + ' [m/s], '                 \
               + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '                   \
               + 'Launch declination = ' + str( launchDeclination[ 0 ] ) + ' [deg]\n',               \
               fontsize=10 )

ax1.grid( )
limits = max( max( x.tolist(), y.tolist(), z.tolist() ) )
ax1.set_xlim( [ -limits, limits ] )
ax1.set_ylim( [ -limits, limits ] )
ax1.set_zlim( [ -0.5e5, 0.5e5 ] )
# ax1.axis('equal')

## plot the 3d trajectory in inertial frame
ax2.plot( inertial_x, inertial_y, inertial_z, zdir = 'z', color=colors.cnames["purple"], linewidth=1.0 )

ax2.text( inertial_x[0], inertial_y[0], inertial_z[0], 'start', size=10, zorder=1, color=colors.cnames["black"] )
ax2.text( inertial_x[endIndex-1], inertial_y[endIndex-1], inertial_z[endIndex-1],
          'end', size=10, zorder=1, color=colors.cnames["black"] )

ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.set_zlabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
# ax2.set_title( 'Particle trajectory around asteroid Eros (Inertial frame) \n'
#                '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '                 \
#                'Launch azimuth=' + str( launchAzimuth[ 0 ] ) + '[deg], '                    \
#                'time=' + str( t[ endIndex-1 ] / (24.0*60.0*60.0) ) + '[day(s)]' )

# ax2.set_title( 'Ellipsoid longest edge, Inertial frame trajectory\n'
#                + '$V_{launch}$ = ' + str( velocityMagnitude[ 0 ] ) + ' [m/s], '                 \
#                + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '                   \
#                + 'Solar phase = ' + str( solarPhase ) + ' [deg]\n'                              \
#                + 'Time = ' + str( lowerTime ) + ' to '                                          \
#                + str( upperTime / (24.0 * 60.0 * 60.0) ) + ' [days]',                           \
#                fontsize=10 )

ax2.set_title( 'Ellipsoid longest edge, Inertial frame trajectory\n'
               + '$V_{launch}$ = ' + str( velocityMagnitude[ 0 ] ) + ' [m/s], '                 \
               + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '                   \
               + 'Launch declination = ' + str( launchDeclination[ 0 ] ) + ' [deg]\n',               \
               fontsize=10 )

ax2.grid( )
limits = max( max( inertial_x.tolist(), inertial_y.tolist(), inertial_z.tolist() ) )
ax2.set_xlim( [ -limits, limits ] )
ax2.set_ylim( [ -limits, limits ] )
ax2.set_zlim( [ -0.5e5, 0.5e5 ] )
# ax2.axis('equal')

## plot the 2D trajectory projections now
# get the end point for a data array
endIndex = np.size( x )

string1 = 'Particle trajectory projection around asteroid Eros (Body fixed frame) \n'
string2 = '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '
string3 = 'Launch azimuth=' + str( launchAzimuth[ 0 ] ) + '[deg], '
string4 = 'time=' + str( t[ endIndex-1 ] / (60.0*60.0) ) + '[hrs]'
topTitleBodyFrame = string1 + string2 + string3 + string4

string1 = 'Particle trajectory projection around asteroid Eros (Inertial frame) \n'
string2 = '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '
string3 = 'Launch azimuth=' + str( launchAzimuth[ 0 ] ) + '[deg], '
string4 = 'time=' + str( t[ endIndex-1 ] / (60.0*60.0) ) + '[hrs]'
topTitleInertialFrame = string1 + string2 + string3 + string4

## Set up the figure
fig = plt.figure( figsize=( 7, 10 ) )
plt.suptitle( 'Ellipsoid longest edge, Rotating frame trajectory\n'
               + '$V_{launch}$ = ' + str( velocityMagnitude[ 0 ] ) + ' [m/s], '                 \
               + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '                   \
               + 'Solar phase = ' + str( solarPhase ) + ' [deg]\n'                              \
               + 'Time = ' + str( lowerTime ) + ' to '                                          \
               + str( upperTime / (24.0 * 60.0 * 60.0) ) + ' [days]',                           \
               fontsize=10 )
gs = gridspec.GridSpec( 3, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

theta = np.linspace(0, 2 * np.pi, 360)

## Common format parameters
ellipsoidColor  = colors.cnames["slategray"]
trajectoryColor = colors.cnames["purple"]
textColor       = colors.cnames["black"]
startColor      = colors.cnames["darkgreen"]
endColor        = colors.cnames["darkred"]

## XY Projection
plotEllipse( alpha, beta, theta, ax1 )

ax1.plot( x, y, color=trajectoryColor, linewidth=0.2 )

## indicate starting point
ax1.text( x[0], y[0], 'start', size=12, color=startColor )

## indicate ending point
ax1.text( x[endIndex-1], y[endIndex-1], 'end', size=12, color=endColor )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.grid( True )
# ax1.axis('equal')

## YZ Projection
plotEllipse( beta, gamma, theta, ax2 )

ax2.plot( y, z, color=trajectoryColor, linewidth=0.2 )

## indicate starting point
ax2.text( y[0], z[0], 'start', size=12, color=startColor, rotation='vertical', va='top' )

## indicate ending point
ax2.text( y[endIndex-1], z[endIndex-1], 'end', size=12, color=endColor, rotation='vertical' )

## format axis and title
ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax2.grid( True )
# ax2.axis('equal')

## XZ Projection
plotEllipse( alpha, gamma, theta, ax3 )

ax3.plot( x, z, color=trajectoryColor, linewidth=0.2 )

## indicate starting point
ax3.text( x[0], z[0], 'start', size=12, color=startColor, rotation='vertical', va='top' )

## indicate ending point
ax3.text( x[endIndex-1], z[endIndex-1], 'end', size=12, color=endColor, rotation='vertical' )

## format axis and title
ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax3.grid( True )
# ax3.axis('equal')


## Inertial Frame
## Set up the figure
fig = plt.figure( figsize=( 7, 10 ) )
plt.suptitle( 'Ellipsoid longest edge, Inertial frame trajectory\n'
               + '$V_{launch}$ = ' + str( velocityMagnitude[ 0 ] ) + ' [m/s], '                 \
               + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '                   \
               + 'Solar phase = ' + str( solarPhase ) + ' [deg]\n'                              \
               + 'Time = ' + str( lowerTime ) + ' to '                                          \
               + str( upperTime / (24.0 * 60.0 * 60.0) ) + ' [days]',                           \
               fontsize=10 )
gs = gridspec.GridSpec( 3, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

## XY Projection
plotEllipse( alpha, beta, theta, ax1 )

ax1.plot( inertial_x, inertial_y, color=trajectoryColor, linewidth=0.2 )

## indicate starting point
ax1.text( inertial_x[0], inertial_y[0], 'start', size=12, color=startColor )

## indicate ending point
ax1.text( inertial_x[endIndex-1], inertial_y[endIndex-1], 'end', size=12, color=endColor )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.grid( True )
# ax1.axis('equal')

## YZ Projection
plotEllipse( beta, gamma, theta, ax2 )

ax2.plot( inertial_y, inertial_z, color=trajectoryColor, linewidth=0.2 )

## indicate starting point
ax2.text( inertial_y[0], inertial_z[0], 'start', size=12, color=startColor, rotation='vertical', va='top' )

## indicate ending point
ax2.text( inertial_y[endIndex-1], inertial_z[endIndex-1], 'end', size=12, color=endColor, rotation='vertical' )

## format axis and title
ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax2.grid( True )
# ax2.axis('equal')

## XZ Projection
plotEllipse( alpha, gamma, theta, ax3 )

ax3.plot( inertial_x, inertial_z, color=trajectoryColor, linewidth=0.2 )

## indicate starting point
ax3.text( inertial_x[0], inertial_z[0], 'start', size=12, color=startColor, rotation='vertical', va='top' )

## indicate ending point
ax3.text( inertial_x[endIndex-1], inertial_z[endIndex-1], 'end', size=12, color=endColor, rotation='vertical' )

## format axis and title
ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax3.grid( True )
# ax3.axis('equal')

## seperate figure for the inertial frame 3D trajectory
## plot the 3d trajectory in inertial frame
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax2 = plt.subplot( gs[ 0 ], projection = '3d' )
ax2.plot( inertial_x, inertial_y, inertial_z, zdir = 'z', color=trajectoryColor, linewidth=0.5 )

ax2.text( inertial_x[0], inertial_y[0], inertial_z[0], 'start', size=10, zorder=1, color='black' )
ax2.text( inertial_x[endIndex-1], inertial_y[endIndex-1], inertial_z[endIndex-1],
          'end', size=10, zorder=1, color='black' )

ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.set_zlabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax2.set_title( 'Ellipsoid longest edge, Inertial frame trajectory\n'
               + '$V_{launch}$ = ' + str( velocityMagnitude[ 0 ] ) + ' [m/s], '                 \
               + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '                   \
               + 'Solar phase = ' + str( solarPhase ) + ' [deg]\n'                              \
               + 'Time = ' + str( lowerTime ) + ' to '                                          \
               + str( upperTime / (24.0 * 60.0 * 60.0) ) + ' [days]',                           \
               fontsize=10 )
ax2.grid( )
# ax2.axis('equal')

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
