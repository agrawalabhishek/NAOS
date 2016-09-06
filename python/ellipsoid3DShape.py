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
# print x[1:100]
y = data['Y'].values
# print y
z = data['Z'].values
# print z
r = np.sqrt( x**2 + y**2 + z**2 )

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

x1 = np.linspace( 0, x.max() )
y1 = np.linspace( 0, y.max() )
xx, yy = np.meshgrid( x1, y1 )
temp = 7000.0**2 * ( 1.0 - ( ( xx**2 / 20000.0**2 ) + ( yy**2 / 7000.0**2 ) ) )
zz = np.sqrt( temp )
ax.plot_surface( xx, yy, zz, rstride=2, cstride=2, cmap=cm.jet, linewidth=0.1, antialiased=False )
plt.show()

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
