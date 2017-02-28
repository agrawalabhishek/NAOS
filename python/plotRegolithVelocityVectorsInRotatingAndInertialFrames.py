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

## Operations
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593
mu = 876514

## Operations
# Connect to SQLite database.
try:
    # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db" )
    database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

phaseAngle = 'N.A.'

fig = plt.figure( )
gs  = gridspec.GridSpec( 1, 2 )
ax1 = plt.subplot( gs[ 0 ], projection = '3d' )
ax2 = plt.subplot( gs[ 1 ], projection = '3d' )

## 3D figure, ellipsoid with velocity vectors at certain angles
launch_angles = range( 0, 360, 90 )
launch_angles.append( 30 )
launch_angle_tuple = tuple( launch_angles )

data1 = pd.read_sql( "SELECT    initial_position_x,                                             \
                                initial_position_y,                                             \
                                initial_position_z,                                             \
                                initial_velocity_x,                                             \
                                initial_velocity_y,                                             \
                                initial_velocity_z,                                             \
                                initial_inertial_velocity_x,                                    \
                                initial_inertial_velocity_y,                                    \
                                initial_inertial_velocity_z,                                    \
                                ROUND( launch_azimuth )                                         \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( launch_azimuth ) IN " + str( launch_angle_tuple ) + "    \
                     AND        start_flag = 1                                                  \
                     AND        ROUND( initial_velocity_magnitude ) = 7;",                      \
                     database )

data1.columns = [ 'pos_x',                                                                  \
                  'pos_y',                                                                  \
                  'pos_z',                                                                  \
                  'vel_x',                                                                  \
                  'vel_y',                                                                  \
                  'vel_z',                                                                  \
                  'inertial_vel_x',                                                         \
                  'inertial_vel_y',                                                         \
                  'inertial_vel_z',                                                         \
                  'launch_azimuth' ]

pos_x                           = data1[ 'pos_x' ]
pos_y                           = data1[ 'pos_y' ]
pos_z                           = data1[ 'pos_z' ]
vel_x                           = data1[ 'vel_x' ]
vel_y                           = data1[ 'vel_y' ]
vel_z                           = data1[ 'vel_z' ]
inertial_vel_x                  = data1[ 'inertial_vel_x' ]
inertial_vel_y                  = data1[ 'inertial_vel_y' ]
inertial_vel_z                  = data1[ 'inertial_vel_z' ]
azimuth                         = data1[ 'launch_azimuth' ]

# plot the ellipsoid
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

newColor = colors.cnames["black"]
surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         rstride=5, cstride=5, alpha=0.3 )
surf.set_facecolor( newColor )
surf.set_linewidth( 0.1 )

surf = ax2.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         rstride=5, cstride=5, alpha=0.3 )
surf.set_facecolor( newColor )
surf.set_linewidth( 0.1 )

# plot the velocity vector now
colors = plt.cm.Vega20( np.linspace( 0, 1, len( azimuth ) ) )

for index in range( 0, len( azimuth ) ):
    ax1.quiver3D( pos_x[ index ], pos_y[ index ], pos_z[ index ],
                  vel_x[ index ], vel_y[ index ], vel_z[ index ],
                  length=500, lw=1, pivot='tail', arrow_length_ratio=0.2,
                  color=colors[ index ], linestyles='solid',
                  label="Azimuth = " + str( azimuth[ index ] ) + " [deg]" )

    ax2.quiver3D( pos_x[ index ], pos_y[ index ], pos_z[ index ],
                  inertial_vel_x[ index ], inertial_vel_y[ index ], inertial_vel_z[ index ],
                  length=500, lw=1, pivot='tail', arrow_length_ratio=0.2,
                  color=colors[ index ], linestyles='solid',
                  label="Azimuth = " + str( azimuth[ index ] ) + " [deg]" )

ax1.legend( ).draggable( )
ax1.grid( )
ax1.set_title( 'Rotating frame' )
ax1.set_xlabel( 'x [m]' )
ax1.set_ylabel( 'y [m]' )
ax1.set_zlabel( 'z [m]' )

ax2.legend( ).draggable( )
ax2.grid( )
ax2.set_title( 'Inertial frame' )
ax2.set_xlabel( 'x [m]' )
ax2.set_ylabel( 'y [m]' )
ax2.set_zlabel( 'z [m]' )

## set global plot title
plt.suptitle( 'Launch velocity vectors for ejecta material \n Longest edge, $V_{mag}$ = 7.0 [m/s], \
              Phase angle = ' + str( phaseAngle ) )

## Plot the metadata (velocity vector components)
fig = plt.figure( )
gs  = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax1.axis( 'off' )

metadata_table = []

columnLabels = [ "Launch \n Azimuth [deg]",             \
                 "Component",                           \
                 "V (Rotating \n Frame) [m/s]",         \
                 "W x R [m/s]",                         \
                 "V (Inertial \n Frame) [m/s]" ]

OmegaCrossPosition_Y = Wz * alpha

# add data for all azimuths
for index in range( 0, len( azimuth ) ):
    metadata_table.append( [ azimuth[ index ], "x", vel_x[ index ], 0.0, inertial_vel_x[ index ] ] )
    metadata_table.append( [ azimuth[ index ], "y", vel_y[ index ], OmegaCrossPosition_Y, inertial_vel_y[ index ] ] )
    metadata_table.append( [ azimuth[ index ], "z", vel_z[ index ], 0.0, inertial_vel_z[ index ] ] )

table = ax1.table( cellText = metadata_table,
                   colLabels = columnLabels,
                   cellLoc = 'center', loc = 'center' )

table.auto_set_font_size(False)
table.set_fontsize( 12 )
table_properties = table.properties( )
table_cells = table_properties[ 'child_artists' ]
for cell in table_cells: cell.set_height( 0.10 )
cell_dict = table.get_celld( )
# for row in xrange( 0, 3 * len( azimuth ) ): cell_dict[ ( row, 2 ) ].set_width( 0.1 )

##close the database
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
