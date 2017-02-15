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
print "--------------------------------------------------------------------------------"
print "                                 NAOS                                           "
print "                                                                                "
print "         Copyright (c) 2016, A. Agrawal (abhishek.agrawal@protonmail.com)       "
print "--------------------------------------------------------------------------------"
print ""


def plotEllipse( semiMajor, semiMinor, angleRange, plotHandle ):
        r = semiMajor * semiMinor / np.sqrt( ( semiMinor * np.cos( angleRange ) )**2
                                           + ( semiMajor * np.sin( angleRange ) )**2 )
        x = r * np.cos( angleRange )
        y = r * np.sin( angleRange )
        plotHandle.plot( x, y )
        plotHandle.grid( )
        plotHandle.hold( True )

## function to convert keplerian elements into inertial coordinates
# All angles in the arguments have to be provided in radians (VERY IMPORTANT)
def convertKeplerElementsToCartesianCoordinates( sma,
                                                 eccentricity,
                                                 inclination,
                                                 raan,
                                                 aop,
                                                 ta,
                                                 mu,
                                                 tolerance ):
    # function defination begins here
    if np.abs( eccentricity - 1.0 ) > tolerance:
        semiLatus = sma * ( 1.0 - eccentricity**2 )
    # for parabolic orbits
    else:
        semiLatus = sma

    # compute position and velocity in the perifocal coordinate system
    radius = semiLatus / ( 1.0 + eccentricity * np.cos( ta ) )

    xPositionPerifocal = radius * np.cos( ta )
    yPositionPerifocal = radius * np.sin( ta )

    xVelocityPerifocal = -1.0 * np.sqrt( mu / semiLatus ) * np.sin( ta )
    yVelocityPerifocal = np.sqrt( mu / semiLatus ) * ( eccentricity + np.cos( ta ) )

    # Calculate components of the rotation matrix
    rotation11 = np.cos( raan ) * np.cos( aop ) - np.sin( raan ) * np.sin( aop ) * np.cos( inclination )

    rotation12 = -1.0 * np.cos( raan ) * np.sin( aop ) - np.sin( raan ) * np.cos( aop ) * np.cos( inclination )

    rotation21 = np.sin( raan ) * np.cos( aop ) + np.cos( raan ) * np.sin( aop ) * np.cos( inclination )

    rotation22 = -1.0 * np.sin( raan ) * np.sin( aop ) + np.cos( raan ) * np.cos( aop ) * np.cos( inclination )

    rotation31 = np.sin( aop ) * np.sin( inclination )

    rotation32 = np.cos( aop ) * np.sin( inclination )

    # Compute the cartesian position and velocity in the inertial frame
    xPosition = rotation11 * xPositionPerifocal + rotation12 * yPositionPerifocal

    yPosition = rotation21 * xPositionPerifocal + rotation22 * yPositionPerifocal

    zPosition = rotation31 * xPositionPerifocal + rotation32 * yPositionPerifocal

    xVelocity = rotation11 * xVelocityPerifocal + rotation12 * yVelocityPerifocal

    yVelocity = rotation21 * xVelocityPerifocal + rotation22 * yVelocityPerifocal

    zVelocity = rotation31 * xVelocityPerifocal + rotation32 * yVelocityPerifocal

    # return the values as a tuple
    return ( xPosition, yPosition, zPosition )

# Start timer.
start_time = time.time( )

# Connect to SQLite database.
try:
        database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

data = pd.read_sql( "SELECT     sma,                                                \
                                eccentricity,                                       \
                                inclination,                                        \
                                raan,                                               \
                                aop,                                                \
                                ta,                                                 \
                                position_x,                                         \
                                position_y,                                         \
                                position_z,                                         \
                                ROUND( initial_velocity_magnitude ),                \
                                inertial_position_x,                                \
                                inertial_position_y,                                \
                                inertial_position_z,                                \
                                ROUND( launch_azimuth ),                            \
                                time                                                \
                     FROM       regolith_trajectory_results                         \
                     WHERE      ROUND( launch_azimuth ) = 120.0                     \
                     AND        ROUND( initial_velocity_magnitude ) = 1             \
                     AND        start_flag = 1;",                                   \
                     database )

data.columns = [ 'sma',                                                 \
                 'eccentricity',                                        \
                 'inclination',                                         \
                 'raan',                                                \
                 'aop',                                                 \
                 'ta',                                                  \
                 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'velocity_magnitude',                                  \
                 'inertial_x',                                          \
                 'inertial_y',                                          \
                 'inertial_z',                                          \
                 'launch_azimuth',                                      \
                 'time' ]

sma                 = data[ 'sma' ]
eccentricity        = data[ 'eccentricity' ]
inclination         = data[ 'inclination' ]
raan                = data[ 'raan' ]
aop                 = data[ 'aop' ]
ta                  = data[ 'ta' ]
x                   = data[ 'x' ]
y                   = data[ 'y' ]
z                   = data[ 'z' ]
velocityMagnitude   = data[ 'velocity_magnitude' ]
inertial_x          = data[ 'inertial_x' ]
inertial_y          = data[ 'inertial_y' ]
inertial_z          = data[ 'inertial_z' ]
launchAzimuth       = data[ 'launch_azimuth' ]
t                   = data[ 'time' ]

inclination = inclination * np.pi / 180.0
raan        = raan * np.pi / 180.0
aop         = aop * np.pi / 180.0

sma = sma.tolist()
eccentricity = eccentricity.tolist()
inclination = inclination.tolist()
raan = raan.tolist()
aop = aop.tolist()

## Set up the figure
fig = plt.figure( )
# get the end point for a data array
endIndex = np.size( x )
plt.suptitle( 'Initial orbital elements trajectory (Inertial frame) \n'
              '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '                 \
              'Launch azimuth=' + str( launchAzimuth[ 0 ] ) + '[deg]' )

gs = gridspec.GridSpec( 3, 1, height_ratios = [ 1, 1, 1 ] )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

# convert the orbital elements into cartesian coordinates for all true anomalies
# ta = range( 0.0, 360.0, 1.0 )
ta = np.linspace( 0.0, 359.0, 1000 )
ta = ta * np.pi / 180.0

tolerance = 1.0e-15

# grav parameter of the asteroid/ellipsoid
mu = 876514

xCoordinate = []
yCoordinate = []
zCoordinate = []

for index in range( 0, len( ta ) ):
    ( xPosition, yPosition, zPosition ) = convertKeplerElementsToCartesianCoordinates(
                                                        sma[0],
                                                        eccentricity[0],
                                                        inclination[0],
                                                        raan[0],
                                                        aop[0],
                                                        ta[index],
                                                        mu,
                                                        tolerance )

    xCoordinate.append( xPosition )
    yCoordinate.append( yPosition )
    zCoordinate.append( zPosition )

## ellipsoidal shape model parameters for the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593

theta = np.linspace(0, 2 * np.pi, 300)

## Common format parameters
ellipsoidColor  = colors.cnames["slategray"]
trajectoryColor = colors.cnames["purple"]
textColor       = colors.cnames["black"]
startColor      = colors.cnames["darkgreen"]
endColor        = colors.cnames["darkred"]

###############################################################
######################## XY Projection ########################

plotEllipse( alpha, beta, theta, ax1 )

ax1.plot( xCoordinate, yCoordinate, color=trajectoryColor )

# ## indicate starting point
# ax1.text( x[0], y[0], 'start', size=12, color=startColor )

# ## indicate ending point
# ax1.text( x[endIndex-1], y[endIndex-1], 'end', size=12, color=endColor )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.grid( True )

###############################################################
######################## YZ Projection ########################

plotEllipse( beta, gamma, theta, ax2 )

ax2.plot( yCoordinate, zCoordinate, color=trajectoryColor )

# ## indicate starting point
# ax2.text( y[0], z[0], 'start', size=12, color=startColor, rotation='vertical', va='top' )

# ## indicate ending point
# ax2.text( y[endIndex-1], z[endIndex-1], 'end', size=12, color=endColor, rotation='vertical' )

## format axis and title
ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax2.grid( True )

###############################################################
######################## XZ Projection ########################

plotEllipse( alpha, gamma, theta, ax3 )

ax3.plot( xCoordinate, zCoordinate, color=trajectoryColor )

# ## indicate starting point
# ax3.text( x[0], z[0], 'start', size=12, color=startColor, rotation='vertical', va='top' )

# ## indicate ending point
# ax3.text( x[endIndex-1], z[endIndex-1], 'end', size=12, color=endColor, rotation='vertical' )

## format axis and title
ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax3.grid( True )

## Show the plot
# plt.tight_layout( )
# plt.grid( True )
plt.show( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
