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
from matplotlib import animation
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

# draw the ellipse using the patches package
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

## static single ellipse drawing (works)
acwRotationAngle = 0.0
ellipse = mpatches.Ellipse( ( 0.0, 0.0 ), 40.0, 14.0, acwRotationAngle, alpha=0.5, fill=False,
                            edgecolor='red' )
ax1.add_patch( ellipse )

def init( ):
    ax1.plot( -20.0, 0.0 )
    # patches = []
    # acwRotationAngle = 0.0
    # ellipse = mpatches.Ellipse( ( 0.0, 0.0 ), 40.0, 14.0, acwRotationAngle, alpha=0.5, fill=False,
    #                             edgecolor='red' )
    # # patches.append( ax1.add_patch( ellipse ) )
    # patches = ( ax1.add_patch( ellipse ) )
    # return patches
    return ellipse

def animate( i ):
    point.set_data( coordinate[ 0 ][ i ], coordinate[ 1 ][ i ] )
    # patches = []
    # acwRotationAngle = RotationAngles[ i ]
    # ellipse = mpatches.Ellipse( ( 0.0, 0.0 ), 40.0, 14.0, acwRotationAngle, alpha=0.5, fill=False,
    #                             edgecolor='red' )
    # # patches.append( ax1.add_patch( ellipse ) )
    # patches = ( ax1.add_patch( ellipse ) )
    # return patches, point
    ellipse.angle = RotationAngles[ i ]
    return ellipse, point

## main segment of the code
# define the rotation angles for the ellipse
RotationAngles = np.arange( 0.0, 360.0, 10.0 )

# define the scatter point and the coordinates for which it has to be animated
point, = ax1.plot( [], [], marker='o' )
y = np.linspace( -20.0, 20.0, len(RotationAngles) )
x = 20.0 * np.ones( len( y ) )
coordinate = [ x, y ]
anim = animation.FuncAnimation( fig, animate, init_func=init,
                                frames=len( RotationAngles ), interval=500, blit=False )

ax1.grid(True)
ax1.set_xlim( -40.0, 40.0 )
ax1.set_ylim( -20.0, 20.0 )
ax1.set_aspect( 1 )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## show the plot
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
