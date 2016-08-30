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

# x = data['X'].values
# y = data['Y'].values
# z = data['Z'].values

# get linearly spaced data between the minimum and maximum values. unique() returns an array of
# unique values in the panda series.
x1 = np.linspace( data['X'].min(), data['X'].max(), len(data['X'].unique( ) ) )
y1 = np.linspace( data['Y'].min(), data['Y'].max(), len(data['Y'].unique( ) ) )
x2, y2 = np.meshgrid( x1, y1 )
z2 = griddata( ( data['X'].values, data['Y'].values ), data['Z'].values, (x2, y2), method='cubic' )

surface = ax.plot_surface( x2, y2, z2, rstride=10, cstride=10, cmap=cm.coolwarm, edgecolor='b' )
fig.colorbar(surface, shrink=0.5, aspect=5)

plt.grid()
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
