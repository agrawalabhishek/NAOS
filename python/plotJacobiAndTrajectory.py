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
# Connect to SQLite database.
try:
        database = sqlite3.connect("../data/regolith_launched_from_leading_edge_test/leadingEdgeTest.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

angleArray = range( 0, 1, 1 )

for angleValue in angleArray:
    data = pd.read_sql( "SELECT     position_x,                                                 \
                                    position_y,                                                 \
                                    position_z,                                                 \
                                    velocity_x,                                                 \
                                    velocity_y,                                                 \
                                    velocity_z,                                                 \
                                    time,                                                       \
                                    jacobi_integral                                             \
                         FROM       regolith_trajectory_results                                 \
                         WHERE      ROUND( launch_azimuth ) = " + str(angleValue) + ";",        \
                         database )

    data.columns = [ 'x',                                                   \
                     'y',                                                   \
                     'z',                                                   \
                     'vx',                                                  \
                     'vy',                                                  \
                     'vz',                                                  \
                     'time',                                                \
                     'jacobi_integral' ]

    x = data[ 'x' ]
    y = data[ 'y' ]
    z = data[ 'z' ]
    vx = data[ 'vx' ]
    vy = data[ 'vy' ]
    vz = data[ 'vz' ]
    t = data[ 'time' ]
    jacobi = data[ 'jacobi_integral' ]

    ## Set up the figure
    fig = plt.figure( )
    gs = gridspec.GridSpec( 3, 1, height_ratios = [ 1.5, 2.5, 1.5 ] )
    ax0 = plt.subplot( gs[ 0 ] )
    ax1 = plt.subplot( gs[ 1 ], projection = '3d' )
    ax2 = plt.subplot( gs[ 2 ] )

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
    # surf.set_facecolor( ( 0, 0, 1, 0.5 ) )
    surf.set_facecolor( newColor )
    surf.set_linewidth( 0.1 )

    ax1.hold( True )

    ## Plot 3D trajectory of the orbiting particle
    ax1.plot( x, y, z, zdir = 'z', color=colors.cnames["purple"] )

    ## indicate starting point
    ax1.text( x[0], y[0], z[0], 'start', size=10, zorder=1, color=colors.cnames["black"] )

    ## indicate ending point
    endIndex = np.size( x )
    ax1.text( x[endIndex-1], y[endIndex-1], z[endIndex-1],
              'end', size=10, zorder=1,
              color=colors.cnames["black"] )

    ## format axis and title
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    ax1.set_zlabel('z [m]')
    ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax1.set_title( 'Particle trajectory around asteroid Eros (Body frame)' )

    ## Plot the metadata (initial state vector)
    ax2.axis( 'off' )
    metadata_table = []
    metadata_table.append( [ "Initial X coordinate", x[0], "[m]" ] )
    metadata_table.append( [ "Initial Y coordinate", y[0], "[m]" ] )
    metadata_table.append( [ "Initial Z coordinate", z[0], "[m]" ] )
    metadata_table.append( [ "Initial X velocity", vx[0], "[m/s]" ] )
    metadata_table.append( [ "Initial Y velocity", vy[0], "[m/s]" ] )
    metadata_table.append( [ "Initial Z velocity", vz[0], "[m/s]" ] )
    metadata_table.append( [ "Simulation time", t[endIndex-1], "[s]" ] )
    table = ax2.table( cellText = metadata_table, colLabels = None, cellLoc = 'center', loc = 'center' )
    table_properties = table.properties( )
    table_cells = table_properties[ 'child_artists' ]
    for cell in table_cells: cell.set_height( 0.15 )
    cell_dict = table.get_celld( )
    for row in xrange( 0, 7 ): cell_dict[ ( row, 2 ) ].set_width( 0.1 )

    ## plot the jacobi integral value
    ax0.plot( t, jacobi, color=colors.cnames['purple'], label= 'launch azimuth=' + str(angleValue) )
    ax0.set_xlabel('time [s]')
    ax0.set_ylabel('Jacobi integral')
    ax0.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
    ax0.grid( )
    ax0.legend( bbox_to_anchor=[0.815, 1.15], loc='center' )

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
