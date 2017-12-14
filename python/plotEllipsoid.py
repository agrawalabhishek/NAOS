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

start_time = time.time( )

## Plot the ellipsoidal shape of the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593

fig = plt.figure( figsize=( 6, 6 ) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ], projection = '3d' )
# ax1.set_aspect('equal')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

newColor = colors.cnames["brown"]
ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         rstride=5, cstride=5, color=newColor, linewidth=10.0, alpha=0.7 )

# draw the body frame x-axis
ax1.quiver3D( 0.0, 0.0, 0.0,
              alpha+8000.0, 0.0, 0.0,
              length=1.0, lw=1, pivot='tail', arrow_length_ratio=0.1,
              color=colors.cnames["black"], linestyles='dashed' )

# draw the body frame y-axis
ax1.quiver3D( 0.0, 0.0, 0.0,
              0.0, beta+8000.0, 0.0,
              length=1.0, lw=1, pivot='tail', arrow_length_ratio=0.1,
              color=colors.cnames["black"], linestyles='dashed' )

# draw the body frame z-axis
ax1.quiver3D( 0.0, 0.0, 0.0,
              0.0, 0.0, gamma+8000.0,
              length=1.0, lw=1, pivot='tail', arrow_length_ratio=0.1,
              color=colors.cnames["black"], linestyles='dashed' )


## format axis and title
ax1.set_xlim( [ -alpha, alpha ] )
ax1.set_ylim( [ -alpha, alpha ] )
ax1.set_zlim( [ -alpha, alpha ] )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)

# ax1.axis('equal')
plt.axis('off')

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
