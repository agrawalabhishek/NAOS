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
import math
from scipy.interpolate import griddata
from numpy import linalg as LA

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
import matplotlib.tri as tri
from matplotlib import rcParams
from matplotlib import cm
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import FancyArrowPatch

print ""
print "---------------------------------------------------------------------------------"
print "                                 NAOS                                            "
print "                                                                                 "
print "         Copyright (c) 2016, A. Agrawal (abhishek.agrawal@protonmail.com)        "
print "---------------------------------------------------------------------------------"
print ""

alpha = 20000.0
beta = 7000.0
gamma = 7000.0

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


def dataExtraction( database ):
    data = pd.read_sql( "SELECT     position_x,                                                   \
                                    position_y,                                                   \
                                    position_z                                                    \
                         FROM       regolith_trajectory_results                                   \
                         WHERE      ROUND( start_flag = 1 )                                       \
                         AND        ROUND( launch_declination ) = 0.0                             \
                         AND        ROUND( initial_velocity_magnitude ) = 6.0;",                  \
                         database )

    data.columns = [ 'x',                                                                         \
                     'y',                                                                         \
                     'z' ]

    if database:
        database.close( )

    x                   = data[ 'x' ]
    y                   = data[ 'y' ]
    z                   = data[ 'z' ]

    return [x, y, z]


def cartesianToLatLong( x, y, z ):
    r = np.sqrt( x**2 + y**2 + z**2 )
    Longitude = np.arctan2( y, x ) * 180.0 / np.pi
    Latitude = np.arcsin( z / r ) * 180.0 / np.pi
    return[Longitude, Latitude]

def surfaceCheck( x, y, z ):
    xTerm = x**2 / alpha**2
    yTerm = y**2 / beta**2
    zTerm = z**2 / gamma**2
    checker = xTerm + yTerm + zTerm - 1.0
    return checker

# Start timer.
start_time = time.time( )

# Connect to SQLite database.
try:
    database1 = sqlite3.connect("../data/VandV/"
                                + "new/"
                                + "launch_location/"
                                + "longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database2 = sqlite3.connect("../data/VandV/"
                                + "new/"
                                + "launch_location/"
                                + "intermediateEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database3 = sqlite3.connect("../data/VandV/"
                                + "new/"
                                + "launch_location/"
                                + "shortestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database4 = sqlite3.connect("../data/VandV/"
                                + "new/"
                                + "launch_location/"
                                + "long30_lat60.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

try:
    database5 = sqlite3.connect("../data/VandV/"
                                + "new/"
                                + "launch_location/"
                                + "long225_lat45.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

[data1_x, data1_y, data1_z] = dataExtraction( database1 )
[data2_x, data2_y, data2_z] = dataExtraction( database2 )
[data3_x, data3_y, data3_z] = dataExtraction( database3 )
[data4_x, data4_y, data4_z] = dataExtraction( database4 )
[data5_x, data5_y, data5_z] = dataExtraction( database5 )

print "Processing data now...\n"

print data1_x[0], data1_y[0], data1_z[0]
print data2_x[0], data2_y[0], data2_z[0]
print data3_x[0], data3_y[0], data3_z[0]
print data4_x[0], data4_y[0], data4_z[0]
print data5_x[0], data5_y[0], data5_z[0]
print "\n"

## print lat long angles from initial state vector
[data1_Long, data1_Lat] = cartesianToLatLong( data1_x[0], data1_y[0], data1_z[0] )
print data1_Long, data1_Lat

[data2_Long, data2_Lat] = cartesianToLatLong( data2_x[0], data2_y[0], data2_z[0] )
print data2_Long, data2_Lat

[data3_Long, data3_Lat] = cartesianToLatLong( data3_x[0], data3_y[0], data3_z[0] )
print data3_Long, data3_Lat

[data4_Long, data4_Lat] = cartesianToLatLong( data4_x[0], data4_y[0], data4_z[0] )
print data4_Long, data4_Lat

[data5_Long, data5_Lat] = cartesianToLatLong( data5_x[0], data5_y[0], data5_z[0] )
print data5_Long, data5_Lat

## Plot the ellipsoidal shape of the asteroid
fig = plt.figure( )
ax1 = fig.gca( projection = '3d' )

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

newColor = colors.cnames["slategray"]
surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         color=newColor,
                         rstride=5, cstride=5, alpha=0.3 )
surf.set_linewidth( 3.0 )

# draw the position vector
positionVector = [ data1_x[0], data1_y[0], data1_z[0] ]
ax1.quiver3D( 0.0, 0.0, 0.0,
              positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              lw = 1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["black"],
              linestyles='solid',
              label='Lat.=$0\degree$, Long.=$0\degree$' )

positionVector = [ data2_x[0], data2_y[0], data2_z[0] ]
ax1.quiver3D( 0.0, 0.0, 0.0,
              positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              lw = 1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["blue"],
              linestyles='solid',
              label='Lat.=$0\degree$, Long.=$90\degree$' )

positionVector = [ data3_x[0], data3_y[0], data3_z[0] ]
ax1.quiver3D( 0.0, 0.0, 0.0,
              positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              lw = 1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["red"],
              linestyles='solid',
              label='Lat.=$90\degree$, Long.=$0\degree$' )

positionVector = [ data4_x[0], data4_y[0], data4_z[0] ]
ax1.quiver3D( 0.0, 0.0, 0.0,
              positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              lw = 1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["orange"],
              linestyles='solid',
              label='Lat.=$60\degree$, Long.=$30\degree$' )

positionVector = [ data5_x[0], data5_y[0], data5_z[0] ]
ax1.quiver3D( 0.0, 0.0, 0.0,
              positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              lw = 1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["purple"],
              linestyles='solid',
              label='Lat.=$45\degree$, Long.=$225\degree$' )

## format axis and title
ax1.set_xlim( [ -20000, 20000 ] )
ax1.set_ylim( [ -20000, 20000 ] )
ax1.set_zlim( [ -20000, 20000 ] )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')
ax1.set_title('Position vector to launch location')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)
ax1.grid(True)
ax1.legend( ).draggable( )
# plt.axis('off')

## check if all launch locations are on the ellipsoid's surface
print "Launch location on CDE surface check...\n"

data1SurfaceCheck = surfaceCheck( data1_x[0], data1_y[0], data1_z[0] )
print data1SurfaceCheck

data2SurfaceCheck = surfaceCheck( data2_x[0], data2_y[0], data2_z[0] )
print data2SurfaceCheck

data3SurfaceCheck = surfaceCheck( data3_x[0], data3_y[0], data3_z[0] )
print data3SurfaceCheck

data4SurfaceCheck = surfaceCheck( data4_x[0], data4_y[0], data4_z[0] )
print data4SurfaceCheck

data5SurfaceCheck = surfaceCheck( data5_x[0], data5_y[0], data5_z[0] )
print data5SurfaceCheck

## Show the plot
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
