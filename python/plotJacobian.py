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
data = pd.read_csv( '../data/regolith_launched_from_leading_edge_test/test.csv' )
x = data['x'].values
y = data['y'].values
z = data['z'].values
vx = data['vx'].values
vy = data['vy'].values
vz = data['vz'].values
t = data['t'].values
jacobian = data['jacobi'].values

# Connect to SQLite database.
# try:
#         database = sqlite3.connect("../data/regolith_launched_from_leading_edge/leadingEdge.db")

# except sqlite3.Error, e:
#         print "Error %s:" % e.args[0]
#         sys.exit(1)

# data = pd.read_sql( "SELECT     jacobi_integral,                        \
#                                 time                                    \
#                      FROM       regolith_trajectory_results             \
#                      WHERE      ROUND( launch_azimuth ) = 100.0;",      \
#                      database )

# data.columns = [ 'jacobi_integral',                                     \
#                  'time' ]

# t = data[ 'time' ]
# jacobian = data[ 'jacobi_integral' ]

# convert time in seconds to earth days
# t = t / ( 24.0 * 60.0 * 60.0 )

# upto 7 significant digits after the decimal point in jacobian
for index in range( 0, len(jacobian) ):
    jacobian[index] = float( "{0:.7f}".format(jacobian[index]) )

## Set up the figure
fig = plt.figure( )
# ax1 = fig.add_subplot( 211, projection = '3d' )
# ax2 = fig.add_subplot( 212, frameon = False )
gs = gridspec.GridSpec( 2, 1, height_ratios = [ 2.5, 1 ] )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )

## plot the jacobian
ax1.plot( t, jacobian, color=colors.cnames['purple'] )
ax1.set_xlabel('time [Earth days]')
ax1.set_ylabel('Jacobi integral')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax1.grid( )

## plot relative change in jacobian
# initialJacobian = jacobian[ 0 ];
# relativeChangeInJacobian = (jacobian - initialJacobian) / initialJacobian

# ax2.plot( t, relativeChangeInJacobian, color=colors.cnames['purple'] )
# ax2.set_xlabel('time [Earth days]')
# ax2.set_ylabel('Relative change in Jacobi integral')
# ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
# ax2.grid( )

## plot meta data
endIndex = np.size( x )
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

## show plot
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
