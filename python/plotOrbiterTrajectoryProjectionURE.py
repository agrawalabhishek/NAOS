'''
Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
Distributed under the MIT License.
See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
'''

# Set up modules and packages
# I/O
import csv
from pprint import pprint

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
data = pd.read_csv( '../data/eomOrbiterURESolution.csv' )
x = data[ 'x' ].values
y = data[ 'y' ].values
z = data[ 'z' ].values
vx = data[ 'vx' ].values
vy = data[ 'vy' ].values
vz = data[ 'vz' ].values
t = data[ 't' ].values

## Set up the figure
fig = plt.figure( )
plt.suptitle( "Particle trajectory projection around asteroid Eros (Body fixed frame)" )
ax1 = fig.add_subplot( 221 )
ax2 = fig.add_subplot( 222 )
ax3 = fig.add_subplot( 223 )
ax4 = fig.add_subplot( 224, frameon=False )

## ellipsoidal shape model parameters for the asteroid
alpha = 20000.0
beta = 20000.0
gamma = 20000.0
Wz = 0.00033118202125129593

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

## Common format parameters
ellipsoidColor = colors.cnames["slategray"]
trajectoryColor = colors.cnames["purple"]
textColor = colors.cnames["black"]
startColor = colors.cnames["darkgreen"]
endColor = colors.cnames["darkred"]

###############################################################
######################## XY Projection ########################

ax1.plot( ellipsoid_x, ellipsoid_y, color=ellipsoidColor )
ax1.hold( True )

ax1.plot( x, y, color=trajectoryColor )

## indicate starting point
ax1.text( x[0], y[0], 'start', size=12, color=startColor )

## indicate ending point
endIndex = np.size( x )
ax1.text( x[endIndex-1], y[endIndex-1], 'end', size=12, color=endColor )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.grid()

###############################################################
######################## YZ Projection ########################

ax2.plot( ellipsoid_y, ellipsoid_z, color=ellipsoidColor )
ax2.hold( True )

ax2.plot( y, z, color=trajectoryColor )

## indicate starting point
ax2.text( y[0], z[0], 'start', size=12, color=startColor, rotation='vertical', va='top' )

## indicate ending point
ax2.text( y[endIndex-1], z[endIndex-1], 'end', size=12, color=endColor, rotation='vertical' )

## format axis and title
ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax2.grid()

###############################################################
######################## XZ Projection ########################

ax3.plot( ellipsoid_x, ellipsoid_z, color=ellipsoidColor )
ax3.hold( True )

ax3.plot( x, z, color=trajectoryColor )

## indicate starting point
ax3.text( x[0], z[0], 'start', size=12, color=startColor, rotation='vertical', va='top' )

## indicate ending point
ax3.text( x[endIndex-1], z[endIndex-1], 'end', size=12, color=endColor, rotation='vertical' )

## format axis and title
ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax3.grid()

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
table_properties = table.properties( )
table_cells = table_properties[ 'child_artists' ]
for cell in table_cells: cell.set_height( 0.10 )
cell_dict = table.get_celld( )
for row in xrange( 0, 7 ): cell_dict[ ( row, 2 ) ].set_width( 0.1 )

####################################################################################################
###################################### Inertial Frame ##############################################

## Set up the figure
fig = plt.figure( )
plt.suptitle( "Particle trajectory projection around asteroid Eros (Inertial frame)" )
ax1 = fig.add_subplot( 221 )
ax2 = fig.add_subplot( 222 )
ax3 = fig.add_subplot( 223 )
ax4 = fig.add_subplot( 224, frameon=False )

## Plot (all of above) with respect to an inertial frame
xInertial = np.zeros( endIndex )
yInertial = np.zeros( endIndex )
zInertial = np.zeros( endIndex )
for i in range(0, endIndex):
    # calculate the angle magnitude first (mod 360)
    phi_Z = Wz * t[i]
    phi_Z = np.mod( phi_Z, 2*np.pi )
    xInertial[i] = x[i]*np.cos(phi_Z) - y[i]*np.sin(phi_Z)
    yInertial[i] = x[i]*np.sin(phi_Z) + y[i]*np.cos(phi_Z)
    zInertial[i] = z[i]

###############################################################
######################## XY Projection ########################

ax1.plot( ellipsoid_x, ellipsoid_y, color=ellipsoidColor )
ax1.hold( True )

ax1.plot( xInertial, yInertial, color=trajectoryColor )

## indicate starting point
ax1.text( xInertial[0], yInertial[0], 'start', size=12, color=startColor )

## indicate ending point
endIndex = np.size( xInertial )
ax1.text( xInertial[endIndex-1], yInertial[endIndex-1], 'end', size=12, color=endColor )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.grid()

###############################################################
######################## YZ Projection ########################

ax2.plot( ellipsoid_y, ellipsoid_z, color=ellipsoidColor )
ax2.hold( True )

ax2.plot( yInertial, zInertial, color=trajectoryColor )

## indicate starting point
ax2.text( yInertial[0], zInertial[0],
          'start', size=12, color=startColor, rotation='vertical', va='top' )

## indicate ending point
ax2.text( yInertial[endIndex-1], zInertial[endIndex-1],
          'end', size=12, color=endColor, rotation='vertical' )

## format axis and title
ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax2.grid()

###############################################################
######################## XZ Projection ########################

ax3.plot( ellipsoid_x, ellipsoid_z, color=ellipsoidColor )
ax3.hold( True )

ax3.plot( xInertial, zInertial, color=trajectoryColor )

## indicate starting point
ax3.text( xInertial[0], zInertial[0],
          'start', size=12, color=startColor, rotation='vertical', va='top' )

## indicate ending point
ax3.text( xInertial[endIndex-1], zInertial[endIndex-1],
          'end', size=12, color=endColor, rotation='vertical' )

## format axis and title
ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax3.grid()

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
table_properties = table.properties( )
table_cells = table_properties[ 'child_artists' ]
for cell in table_cells: cell.set_height( 0.10 )
cell_dict = table.get_celld( )
for row in xrange( 0, 7 ): cell_dict[ ( row, 2 ) ].set_width( 0.1 )

## Show the plot
plt.tight_layout( )
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
