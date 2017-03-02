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
# Read data in csv file. data returned as a panda series.
## data in body frame
# data = pd.read_csv( '../data/trajectory_for_different_launch_azimuth/regolithTrajectoryAtAzimuth210.csv' )
# x = data[ 'x' ].values
# y = data[ 'y' ].values
# z = data[ 'z' ].values
# vx = data[ 'vx' ].values
# vy = data[ 'vy' ].values
# vz = data[ 'vz' ].values
# t = data[ 't' ].values

## Plot the ellipsoidal shape of the asteroid
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
                     AND        ROUND( initial_velocity_magnitude ) = 5.0;",        \
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

## Set up the figure
fig = plt.figure( )
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

newColor = colors.cnames["slategray"]
surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         rstride=5, cstride=5, alpha=0.5 )
# surf.set_facecolor( ( 0, 0, 1, 0.5 ) )
surf.set_facecolor( newColor )
surf.set_linewidth( 0.1 )
ax1.set_xlim( -20000, 20000 )
ax1.set_ylim( -7000, 7000 )
ax1.set_zlim( -7000, 7000 )
ax1.hold( True )

## Plot 3D trajectory of the orbiting particle
ax1.plot( x, y, z, zdir = 'z', color=colors.cnames["purple"] )

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
ax1.set_title( 'Particle trajectory around asteroid Eros (Body frame) \n'
               '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '                 \
               'Launch azimuth=' + str( launchAzimuth[ 0 ] ) + '[deg], '                    \
               'time=' + str( t[ endIndex-1 ] / (24.0*60.0*60.0) ) + '[day(s)]' )
ax1.grid( )
ax1.axis('equal')

## plot the 3d trajectory in inertial frame
ax2.plot( inertial_x, inertial_y, inertial_z, zdir = 'z', color=colors.cnames["purple"] )

ax2.text( inertial_x[0], inertial_y[0], inertial_z[0], 'start', size=10, zorder=1, color=colors.cnames["black"] )
ax2.text( inertial_x[endIndex-1], inertial_y[endIndex-1], inertial_z[endIndex-1],
          'end', size=10, zorder=1, color=colors.cnames["black"] )

ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.set_zlabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax2.set_title( 'Particle trajectory around asteroid Eros (Inertial frame) \n'
               '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '                 \
               'Launch azimuth=' + str( launchAzimuth[ 0 ] ) + '[deg], '                    \
               'time=' + str( t[ endIndex-1 ] / (24.0*60.0*60.0) ) + '[day(s)]' )
ax2.grid( )
ax2.axis('equal')

## Plot the metadata (initial state vector)
# ax2.axis( 'off' )
# metadata_table = []
# metadata_table.append( [ "Initial X coordinate", x[0], "[m]" ] )
# metadata_table.append( [ "Initial Y coordinate", y[0], "[m]" ] )
# metadata_table.append( [ "Initial Z coordinate", z[0], "[m]" ] )
# metadata_table.append( [ "Initial X velocity", vx[0], "[m/s]" ] )
# metadata_table.append( [ "Initial Y velocity", vy[0], "[m/s]" ] )
# metadata_table.append( [ "Initial Z velocity", vz[0], "[m/s]" ] )
# metadata_table.append( [ "Simulation time", t[endIndex-1], "[s]" ] )
# table = ax2.table( cellText = metadata_table, colLabels = None, cellLoc = 'center', loc = 'center' )
# table_properties = table.properties( )
# table_cells = table_properties[ 'child_artists' ]
# for cell in table_cells: cell.set_height( 0.15 )
# cell_dict = table.get_celld( )
# for row in xrange( 0, 7 ): cell_dict[ ( row, 2 ) ].set_width( 0.1 )

## sanity check for the initial position coordinate
ellipsoidSolution = ( np.square(x[0]) / alpha**2 ) + ( np.square(y[0]) / beta**2 ) + ( np.square(z[0]) / gamma**2 )
ellipsoidSolution = ellipsoidSolution - 1.0
if ellipsoidSolution > 1.0e-12:
    print "launch point not on ellipsoid"

## Show the plot
# plt.tight_layout( )
# plt.grid( )
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
