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
        plotHandle.plot( x, y, label='Asteroid at regolith launch' )

# Start timer.
start_time = time.time( )

## Operations
# Read data in csv file. data returned as a panda series.
# data = pd.read_csv( '../data/singleRegolithEjectaURESolution.csv' )
# x = data[ 'x' ].values
# y = data[ 'y' ].values
# z = data[ 'z' ].values
# vx = data[ 'vx' ].values
# vy = data[ 'vy' ].values
# vz = data[ 'vz' ].values
# t = data[ 't' ].values

## ellipsoidal shape model parameters for the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593

# Connect to SQLite database.
try:
        # database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

data = pd.read_sql( "SELECT     position_x,                                         \
                                position_y,                                         \
                                position_z,                                         \
                                ROUND( initial_velocity_magnitude ),                \
                                inertial_position_x,                                \
                                inertial_position_y,                                \
                                inertial_position_z,                                \
                                ROUND( launch_azimuth ),                            \
                                time                                                \
                     FROM       regolith_trajectory_results                         \
                     WHERE      ROUND( launch_azimuth ) = 270.0                     \
                     AND        ROUND( initial_velocity_magnitude ) = 5;",          \
                     database )

data.columns = [ 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'velocity_magnitude',                                  \
                 'inertial_x',                                          \
                 'inertial_y',                                          \
                 'inertial_z',                                          \
                 'launch_azimuth',                                      \
                 'time' ]

x                   = data[ 'x' ]
y                   = data[ 'y' ]
z                   = data[ 'z' ]
velocityMagnitude   = data[ 'velocity_magnitude' ]
inertial_x          = data[ 'inertial_x' ]
inertial_y          = data[ 'inertial_y' ]
inertial_z          = data[ 'inertial_z' ]
launchAzimuth       = data[ 'launch_azimuth' ]
t                   = data[ 'time' ]

t = t.tolist( )

# get the end point for a data array
endIndex = np.size( x )

string1 = 'Particle trajectory (Rotating frame) \n'
string2 = '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '
string3 = 'Launch azimuth=' + str( launchAzimuth[ 0 ] ) + '[deg], '
string4 = 'time=' + str( t[ endIndex-1 ] / (60.0*60.0) ) + '[hrs] \n'
# string5 = 'Phase angle = ' + str( phaseAngle ) + ' [deg]'
topTitleBodyFrame = string1 + string2 + string3 + string4

string1 = 'Particle trajectory (Inertial frame) \n'
string2 = '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '
string3 = 'Launch azimuth=' + str( launchAzimuth[ 0 ] ) + '[deg], '
string4 = 'time=' + str( t[ endIndex-1 ] / (60.0*60.0) ) + '[hrs] \n'
# string5 = 'Phase angle = ' + str( phaseAngle ) + ' [deg]'
topTitleInertialFrame = string1 + string2 + string3 + string4

## Set up the figure
fig = plt.figure( )
plt.suptitle( topTitleBodyFrame )

gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

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

ax1.plot( x, y, color=trajectoryColor )

## indicate starting point
ax1.text( x[0], y[0], 'start', size=12, color=startColor )

## indicate ending point
ax1.text( x[endIndex-1], y[endIndex-1], 'end', size=12, color=endColor )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.grid( True )
ax1.axis('equal')

######################################### Inertial Frame ########################################
#################################################################################################
## Set up the figure
fig = plt.figure( )
plt.suptitle( topTitleInertialFrame )

gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

acwRotationAngle = Wz * t[ -1 ] * 180.0 / np.pi
ellipse2 = mpatches.Ellipse( ( 0.0, 0.0 ), 40.0e3, 14.0e3, acwRotationAngle, alpha=1.0,
                             fill=False,
                             edgecolor='black',
                             label='Asteroid at simulation end' )
ax1.add_patch( ellipse2 )
ax1.plot( [], [] )

ax1.set_aspect( 1 )

###############################################################
######################## XY Projection ########################

plotEllipse( alpha, beta, theta, ax1 )

ax1.plot( inertial_x, inertial_y, color=trajectoryColor, label='Regolith trajectory' )

## indicate starting point
ax1.text( inertial_x[0], inertial_y[0], 'start', size=12, color=startColor )

## indicate ending point
ax1.text( inertial_x[endIndex-1], inertial_y[endIndex-1], 'end', size=12, color=endColor )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.grid( True )
ax1.legend( ).draggable( )
ax1.axis('equal')

## close the database
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
