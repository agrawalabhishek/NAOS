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
from matplotlib.patches import Ellipse

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
fig = plt.figure()
ax1 = fig.add_subplot( 221 )
ax2 = fig.add_subplot( 222 )
ax3 = fig.add_subplot( 223 )

# Read data in csv file. data returned as a panda series.
data = pd.read_csv( '../data/singleRegolithEjectaURESolution.csv' )
x = data[ 'x' ].values
y = data[ 'y' ].values
z = data[ 'z' ].values
vx = data[ 'vx' ].values
vy = data[ 'vy' ].values
vz = data[ 'vz' ].values
t = data[ 't' ].values

## Plot the ellipsoidal projection of the asteroid
alpha = 20000.000
beta = 7000.000
gamma = 7000.000
theta = np.linspace(0, 2 * np.pi, 500)

def plotEllipse( semiMajor, semiMinor, angleRange, plotHandle ):
    r = semiMajor * semiMinor / np.sqrt( ( semiMinor * np.cos( angleRange ) )**2
                                       + ( semiMajor * np.sin( angleRange ) )**2 )
    x = r * np.cos( angleRange )
    y = r * np.sin( angleRange )
    plotHandle.plot( x, y )
    plotHandle.grid( )
    plotHandle.hold( True )

plotEllipse( alpha, beta, theta, ax1 )
plotEllipse( beta, gamma, theta, ax2 )
plotEllipse( alpha, gamma, theta, ax3 )

## draw the vector
velocityVector = [ vx[0], vy[0], vz[0] ]
positionVector = [ x[0], y[0], z[0] ]

## xy plane
# ax1.scatter( positionVector[ 0 ], positionVector[ 1 ], color=colors.cnames["darkred"] )
ax1.quiver( positionVector[ 0 ], positionVector[ 1 ], velocityVector[ 0 ], velocityVector[ 1 ] )

## yz plane
# ax2.scatter( positionVector[ 1 ], positionVector[ 2 ], color=colors.cnames["darkred"] )
ax2.quiver( positionVector[ 1 ], positionVector[ 2 ], velocityVector[ 1 ], velocityVector[ 2 ] )

## xz plane
# ax3.scatter( positionVector[ 0 ], positionVector[ 2 ], color=colors.cnames["darkred"] )
ax3.quiver( positionVector[ 0 ], positionVector[ 2 ], velocityVector[ 0 ], velocityVector[ 2 ] )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)

ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)

ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)

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
