'''
Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
Distributed under the MIT License.
See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
'''

# Set up modules and packages
# I/O
import csv
from pprint import pprint

# 3D visualization special package
import mayavi
from mayavi import mlab

# System
import sys
import time
from tqdm import tqdm

print ""
print "---------------------------------------------------------------------------------"
print "                                 NAOS                                            "
print "                                                                                 "
print "         Copyright (c) 2016, A. Agrawal (abhishek.agrawal@protonmail.com)        "
print "---------------------------------------------------------------------------------"
print ""

# Start timer.
start_time = time.time( )

# Get plotting packages
import matplotlib
import matplotlib.colors
import matplotlib.axes
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import rcParams
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.tri as tri

# Numerical
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import math

# Operations
# Read data in csv file. data returned as a panda series.
data = pd.read_csv( '../data/ellipsoidSurfacePoints.csv' )

# Plot 3D surface of the ellipsoid
fig = plt.figure()
ax = fig.gca( projection = '3d' )
ax.set_xlabel('x [km]')
ax.set_ylabel('y [km]')
ax.set_zlabel('z [km]')
ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

x = data['X'].values
y = data['Y'].values
z = data['Z'].values
latitude = data['latitude'].values
longitude = data['longitude'].values

# **************** trisurf, scatter and wireframe ************ #
# triang = tri.Triangulation( x, y )
# ax.plot_trisurf( x, y, z, triangles=triang.triangles, cmap=cm.jet, linewidth=0.1 )
# ax.scatter( x, y, z )
# ax.plot_wireframe( x, y, z )
# plt.show()

# **************** Method - 1 ******************************** #
# pts = mayavi.mlab.points3d( x, y, z, z )
# mesh = mayavi.mlab.pipeline.delaunay2d( pts )
# pts.remove( )
# surf = mayavi.mlab.pipeline.surface( mesh )
# mayavi.mlab.show( )

# **************** Method - 2 ******************************** #
# x1 = np.linspace( x.min(), x.max() )
# y1 = np.linspace( y.min(), y.max() )
# xx, yy = np.meshgrid( x1, y1 )
# zz = griddata( ( x, y ), z, ( x1, y1 ), method='cubic' )
# ax.plot_surface( xx, yy, zz, rstride=5, cstride=5, cmap=cm.jet, linewidth=0.1, antialiased=False )
# plt.show()

# **************** Method - 3 ******************************** #
# x1 = np.linspace( x.min(), x.max() )
# y1 = np.linspace( y.min(), y.max() )
# xx, yy = np.meshgrid( x1, y1 )
# zz = griddata( ( x, y ), z, ( xx, yy ), method='cubic' )
# ax.plot_surface( xx, yy, zz, rstride=5, cstride=5, cmap=cm.jet, linewidth=0.1, antialiased=False )
# plt.show()

# **************** Method - 4 ******************************** #
# x1 = np.linspace( x.min(), x.max(), len( data['X'].unique() ) )
# y1 = np.linspace( y.min(), y.max(), len( data['Y'].unique() ) )
# xx, yy = np.meshgrid( x1, y1 )
# zz = griddata( ( x, y ), z, ( xx, yy ), method='cubic' )
# ax.plot_surface( xx, yy, zz, rstride=5, cstride=5, cmap=cm.jet, linewidth=0.1, antialiased=False )
# plt.show()

# **************** Method - 5 ******************************** #
# xx, yy = np.mgrid[ min(x):max(x):100j, min(y):max(y):100j ]
# zz = griddata( ( x, y ), z, ( xx, yy ), method='cubic' )
# ax.plot_surface( xx, yy, zz, rstride=5, cstride=5, cmap=cm.jet, linewidth=0.1, antialiased=False )
# plt.show()

# **************** Method - 6 ******************************** #
# x1 = np.linspace( x.min(), x.max() )
# y1 = np.linspace( y.min(), y.max() )
# xx, yy = np.meshgrid( x1, y1 )
# temp = 7000.0**2 * ( 1.0 - ( ( xx**2 / 20000.0**2 ) + ( yy**2 / 7000.0**2 ) ) )
# zz = np.sqrt( temp )
# ax.plot_surface( xx, yy, zz, rstride=2, cstride=2, cmap=cm.jet, linewidth=0.1, antialiased=False )
# plt.show()

# **************** Method - 7 ******************************** #
# cols = np.unique( x ).shape[ 0 ]
# xx = x.reshape( -1, cols )
# yy = y.reshape( -1, cols )
# zz = z.reshape( -1, cols )
# ax.plot_surface( xx, yy, zz, rstride=2, cstride=2, cmap=cm.jet, linewidth=0.1, antialiased=False )
# plt.show()

# **************** Method - 8 ******************************** #
# find number of unique data points
# numberOfLatitudes = len( data['latitude'].unique() )
# numberOfLongitudes = len( data['longitude'].unique() )
# make 2D arrays without changing data
# yy = y.reshape( numberOfLatitudes, numberOfLongitudes )
# xx = x.reshape( numberOfLatitudes, numberOfLongitudes )
# zz = z.reshape( numberOfLatitudes, numberOfLongitudes )
# ax.plot_surface( xx, yy, zz )
# plt.show()

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
