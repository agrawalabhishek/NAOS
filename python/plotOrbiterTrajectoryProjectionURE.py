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

# Connect to SQLite database.
try:
        database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/leadingEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

data = pd.read_sql( "SELECT     position_x,                                         \
                                position_y,                                         \
                                position_z,                                         \
                                velocity_x,                                         \
                                velocity_y,                                         \
                                velocity_z,                                         \
                                time                                                \
                     FROM       regolith_trajectory_results                         \
                     WHERE      ROUND( launch_azimuth ) = 125.0                     \
                     AND        ROUND( initial_velocity_magnitude ) = 14;",         \
                     database )

data.columns = [ 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'vx',                                                  \
                 'vy',                                                  \
                 'vz',                                                  \
                 'time' ]

x = data[ 'x' ]
y = data[ 'y' ]
z = data[ 'z' ]
vx = data[ 'vx' ]
vy = data[ 'vy' ]
vz = data[ 'vz' ]
t = data[ 'time' ]

# mismatchTime = []
# mismatchFlag = False
# for index in range( 1, len(jacobi) ):
#     if np.absolute( jacobi[ index ] - jacobi[ index - 1 ] ) >= 1.0e-10:
#         mismatchFlag = True
#         mismatchTime.append( index )
#     # upto 7 significant digits after the decimal point in jacobian
#     # jacobi[index] = float( "{0:.7f}".format(jacobi[index]) )

## Set up the figure
fig = plt.figure( )
plt.suptitle( "Particle trajectory projection around asteroid Eros (Body fixed frame)" )
ax1 = fig.add_subplot( 221 )
ax2 = fig.add_subplot( 222 )
ax3 = fig.add_subplot( 223 )
ax4 = fig.add_subplot( 224, frameon=False )

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

ax1.plot( x, y, color=trajectoryColor )

## indicate starting point
ax1.text( x[0], y[0], 'start', size=12, color=startColor )

## indicate ending point
endIndex = np.size( x )
ax1.text( x[endIndex-1], y[endIndex-1], 'end', size=12, color=endColor )

# ## plot locations where jacobian switched from its constant value
# ax1.scatter( x[mismatchTime], y[mismatchTime], color='red' )
# for index in range( 0, len( mismatchTime ) ):
#     ax1.text( x[mismatchTime[index]], y[mismatchTime[index]], str( mismatchTime[index] ), color='red' )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

###############################################################
######################## YZ Projection ########################

plotEllipse( beta, gamma, theta, ax2 )

ax2.plot( y, z, color=trajectoryColor )

## indicate starting point
ax2.text( y[0], z[0], 'start', size=12, color=startColor, rotation='vertical', va='top' )

## indicate ending point
ax2.text( y[endIndex-1], z[endIndex-1], 'end', size=12, color=endColor, rotation='vertical' )

# ## plot locations where jacobian switched from its constant value
# ax2.scatter( y[mismatchTime], z[mismatchTime], color='red' )
# for index in range( 0, len( mismatchTime ) ):
#     ax2.text( y[mismatchTime[index]], z[mismatchTime[index]], str( mismatchTime[index] ), color='red' )

## format axis and title
ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

###############################################################
######################## XZ Projection ########################

plotEllipse( alpha, gamma, theta, ax3 )

ax3.plot( x, z, color=trajectoryColor )

## indicate starting point
ax3.text( x[0], z[0], 'start', size=12, color=startColor, rotation='vertical', va='top' )

## indicate ending point
ax3.text( x[endIndex-1], z[endIndex-1], 'end', size=12, color=endColor, rotation='vertical' )

# ## plot locations where jacobian switched from its constant value
# ax3.scatter( x[mismatchTime], z[mismatchTime], color='red' )
# for index in range( 0, len( mismatchTime ) ):
#     ax3.text( x[mismatchTime[index]], z[mismatchTime[index]], str( mismatchTime[index] ), color='red' )

## format axis and title
ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

###############################################################
######################## MetaData #############################

ax4.axis( 'off' )
metadata_table = []
metadata_table.append( [ "Initial X coordinate", x[0], "[m]" ] )
metadata_table.append( [ "Initial Y coordinate", y[0], "[m]" ] )
metadata_table.append( [ "Initial Z coordinate", z[0], "[m]" ] )
metadata_table.append( [ "Initial X velocity", vx[0], "[m/s]" ] )
metadata_table.append( [ "Initial Y velocity", vy[0], "[m/s]" ] )
metadata_table.append( [ "Initial Z velocity", vz[0], "[m/s]" ] )
metadata_table.append( [ "Simulation time", t[endIndex-1], "[s]" ] )
table = ax4.table( cellText = metadata_table, colLabels = None, cellLoc = 'center', loc = 'center' )
table.auto_set_font_size(False)
table.set_fontsize( 12 )
table_properties = table.properties( )
table_cells = table_properties[ 'child_artists' ]
for cell in table_cells: cell.set_height( 0.15 )
cell_dict = table.get_celld( )
for row in xrange( 0, 7 ): cell_dict[ ( row, 2 ) ].set_width( 0.1 )

## Show the plot
plt.tight_layout( )
plt.grid( )
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
