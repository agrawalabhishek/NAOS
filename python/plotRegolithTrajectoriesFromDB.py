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
                         rstride=5, cstride=5 )
surf.set_facecolor( newColor )
surf.set_linewidth( 0.1 )

ax1.hold( True )

## data in body frame
databasePath = "../../data/position_and_launch_angles_iterator/regolithMonteCarlo.db"
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
                     FROM       regolith_trajectory_results             \
                     WHERE      ROUND( launch_azimuth ) = 349.0;",      \
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

## Plot 3D trajectory of the orbiting particle
ax1.plot( x, y, z, zdir = 'z', color=plotColor )

## indicate starting point
ax1.text( x[0], y[0], z[0], 'start', size=10, zorder=1, color=plotColor )

## indicate ending point
endIndex = np.size( x )
ax1.text( x[endIndex-1], y[endIndex-1], z[endIndex-1],
          'end', size=10, zorder=1,
          color=plotColor )

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
