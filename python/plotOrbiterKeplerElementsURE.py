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
# data = pd.read_csv( '../data/eomOrbiterURESolution.csv' )
data = pd.read_csv( '../data/pointMassSolution.csv' )
x = data[ 'x' ].values
y = data[ 'y' ].values
z = data[ 'z' ].values
vx = data[ 'vx' ].values
vy = data[ 'vy' ].values
vz = data[ 'vz' ].values
t = data[ 't' ].values
jacobian = data[ 'jacobian' ].values
semiMajor = data['semiMajor'].values
eccentricity = data['eccentricity'].values
inclination = data['inclination'].values
raan = data['RAAN'].values
aop = data['AOP'].values
ta = data['TA'].values
stepSize = data['stepSize'].values

## Set up the figure
fig = plt.figure( )
ax1 = fig.add_subplot( 421 )
ax2 = fig.add_subplot( 422 )
ax3 = fig.add_subplot( 423 )
ax4 = fig.add_subplot( 424 )
ax5 = fig.add_subplot( 425 )
ax6 = fig.add_subplot( 426 )
ax7 = fig.add_subplot( 427 )

## plot the orbital elements
ax1.plot( t, semiMajor, color=colors.cnames['purple'] )
ax1.set_xlabel('time [s]')
ax1.set_ylabel('semi-major axis [m]')
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
ax1.grid( )

ax2.plot( t, eccentricity, color=colors.cnames['purple'] )
ax2.set_xlabel('time [s]')
ax2.set_ylabel('eccentricity')
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
ax2.grid( )

ax3.plot( t, inclination, color=colors.cnames['purple'] )
ax3.set_xlabel('time [s]')
ax3.set_ylabel('inclination [deg]')
ax3.ticklabel_format(style='plain', useOffset=False)
ax3.grid( )

ax4.plot( t, raan, color=colors.cnames['purple'] )
ax4.set_xlabel('time [s]')
ax4.set_ylabel('RAAN [deg]')
ax4.ticklabel_format(style='plain', useOffset=False)
ax4.grid( )

ax5.plot( t, aop, color=colors.cnames['purple'] )
ax5.set_xlabel('time [s]')
ax5.set_ylabel('AOP [deg]')
ax5.ticklabel_format(style='plain', useOffset=False)
ax5.grid( )

ax6.plot( t, stepSize, color=colors.cnames['purple'] )
ax6.set_xlabel('time [s]')
ax6.set_ylabel('stepSize [s]')
maxStep = np.max( stepSize )
minStep = np.min( stepSize )
ax6.set_ylim( [ -0.01, 0.01 ] )
ax6.ticklabel_format(style='plain', useOffset=False)
ax6.grid( )

## plot meta data
endIndex = np.size( x )
ax7.axis( 'off' )
metadata_table = []
metadata_table.append( [ "Initial X coordinate", x[0], "[m]" ] )
metadata_table.append( [ "Initial Y coordinate", y[0], "[m]" ] )
metadata_table.append( [ "Initial Z coordinate", z[0], "[m]" ] )
metadata_table.append( [ "Initial X velocity", vx[0], "[m/s]" ] )
metadata_table.append( [ "Initial Y velocity", vy[0], "[m/s]" ] )
metadata_table.append( [ "Initial Z velocity", vz[0], "[m/s]" ] )
metadata_table.append( [ "Simulation time", t[endIndex-1], "[s]" ] )
table = ax7.table( cellText = metadata_table, colLabels = None, cellLoc = 'center', loc = 'center' )
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
