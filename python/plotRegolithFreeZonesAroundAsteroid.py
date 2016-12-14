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

## Set up the figure
fig = plt.figure( )
ax1 = fig.gca( projection = '3d' )

## Plot the ellipsoidal shape of the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

newColor = colors.cnames["slategray"]
surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         rstride=5, cstride=5, alpha = 0.5 )
surf.set_facecolor( newColor )
surf.set_linewidth( 0.1 )

ax1.hold( True )

## data in body frame
print "Fetching scan data from database ..."

# Connect to SQLite database.
try:
        database = sqlite3.connect("../data/regolith_launched_from_leading_edge/leadingEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

# data = pd.read_sql( "SELECT     position_x,                             \
#                                 position_y,                             \
#                                 position_z,                             \
#                                 velocity_x,                             \
#                                 velocity_y,                             \
#                                 velocity_z,                             \
#                                 time                                    \
#                      FROM       regolith_trajectory_results             \
#                      WHERE      ROUND( launch_azimuth ) > 208.0         \
#                      AND        ROUND( launch_azimuth ) < 210.0;",      \
#                      database )

data = pd.read_sql( "SELECT     position_x,                             \
                                position_y,                             \
                                position_z,                             \
                                velocity_x,                             \
                                velocity_y,                             \
                                velocity_z,                             \
                                time                                    \
                     FROM       regolith_trajectory_results;",          \
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

plotColor = 'r'

## seperate out points which lie on regolith trajectory within the outer radius
# outerRadius = alpha + 10.0

altitude = 1000.0
# u = np.linspace(0, 2 * np.pi, 200)
# v = np.linspace(0, np.pi, 200)
# outerEllipsoid_x = (alpha+altitude) * np.cos(u) * np.sin(v)
# outerEllipsoid_y = (beta+altitude) * np.sin(u) * np.sin(v)
# outerEllipsoid_z = (gamma+altitude) * np.cos(v)

# outerRadius = []
# outerRadius = np.sqrt( np.ndarray.flatten( outerEllipsoid_x )**2
#     + np.ndarray.flatten( outerEllipsoid_y )**2
#     + np.ndarray.flatten( outerEllipsoid_z )**2 )

# print outerRadius

xNew = []
yNew = []
zNew = []
for index in tqdm( range( len( x ) ) ):
    # radialDistance = np.sqrt( x[index]**2 + y[index]**2 + z[index]**2 )
    xChecker = x[index]**2 / (alpha+altitude)**2
    yChecker = y[index]**2 / (beta+altitude)**2
    zChecker = z[index]**2 / (gamma+altitude)**2
    checker = xChecker + yChecker + zChecker
    if checker <= 1.0:
        xNew.append( x[index] )
        yNew.append( y[index] )
        zNew.append( z[index] )

xSafe = []
ySafe = []
zSafe = []
## get points, within the outerRadius, that are not equal to the points seperated out above
# xIterator_list = np.arange( 0.0, outerRadius+1.0, 5.0 )
# yIterator_list = np.arange( 0.0, outerRadius+1.0, 5.0 )
# zIterator_list = np.arange( 0.0, outerRadius+1.0, 5.0 )
# for xIterator in tqdm( xIterator_list ):
#     for yIterator in yIterator_list:
#         for zIterator in zIterator_list:
#             surfaceCheck = ( xIterator**2 / alpha**2 ) + ( yIterator**2 / beta**2 ) + ( zIterator**2 / gamma**2 )
#             if surfaceCheck >= 1.0: # above or at the surface
#                 radialDistance = np.sqrt( xIterator**2 + yIterator**2 + zIterator**2 )
#                 if radialDistance <= outerRadius:
#                     xFound = xIterator in xNew
#                     yFound = yIterator in yNew
#                     zFound = zIterator in zNew
#                     if xFound and yFound and zFound: # only when all are True, the point is not safe
#                         continue
#                     else:
#                         xSafe.append( xIterator )
#                         ySafe.append( yIterator )
#                         zSafe.append( zIterator )

# u = np.linspace(0, 2 * np.pi, 200)
# v = np.linspace(0, np.pi, 200)

# sma1_list = np.arange( alpha + 1.0, alpha + 11.0, 1.0 )
# sma2_list = np.arange( beta + 1.0, beta + 11.0, 1.0 )
# sma3_list = np.arange( gamma + 1.0, gamma + 11.0, 1.0 )

# for sma1 in tqdm(sma1_list):
#     for sma2 in sma2_list:
#         for sma3 in sma3_list:
#             ellipsoid_x = sma1 * np.outer(np.cos(u), np.sin(v))
#             ellipsoid_y = sma2 * np.outer(np.sin(u), np.sin(v))
#             ellipsoid_z = sma3 * np.outer(np.ones(np.size(u)), np.cos(v))

#             xNew = np.array(list(xNew))
#             ellipsoid_x = np.array(list(ellipsoid_x))
#             ellipsoid_x = np.ndarray.flatten(ellipsoid_x)
#             xSet = set(ellipsoid_x).difference(xNew)
#             xSafe = np.array(list(xSet))

#             yNew = np.array(list(yNew))
#             ellipsoid_y = np.array(list(ellipsoid_y))
#             ellipsoid_y = np.ndarray.flatten(ellipsoid_y)
#             ySet = set(ellipsoid_y).difference(yNew)
#             ySafe = np.array(list(ySet))

#             zNew = np.array(list(zNew))
#             ellipsoid_z = np.array(list(ellipsoid_z))
#             ellipsoid_z = np.ndarray.flatten(ellipsoid_z)
#             zSet = set(ellipsoid_z).difference(zNew)
#             zSafe = np.array(list(zSet))


## Plot all safe points in 3D around the ellipsoid
ax1.scatter( xNew, yNew, zNew, zdir = 'z', color=plotColor )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.set_title( 'Particle trajectory around asteroid Eros (Body frame)' )

## Show the plot
# plt.legend( )
plt.tight_layout( )
plt.grid( )
plt.show( )
# plt.savefig('../data/trajectory_for_different_launch_azimuth/trajectoryPlot_'+str(figureIndex+1)+'.png', dpi=300)

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
