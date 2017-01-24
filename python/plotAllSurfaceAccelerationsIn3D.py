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
from scipy.interpolate import RectBivariateSpline
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

def fmt(x, pos):
    a, b = '{:.7e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

## Operations
# Read data in csv file. data returned as a panda series.
data = pd.read_csv( '../data/surface_gravity_and_perturbing_acceleration/zero_inclination/allAccelerations.csv' )
Ux = data['xGravAcceleration'].values
Uy = data['yGravAcceleration'].values
Uz = data['zGravAcceleration'].values
U = data['gravAcceleration'].values

sunPeriapsisTBP = data['thirdBodyAcceleration_sunPeriapsis'].values
sunPeriapsisSRP = data['srpAcceleration_sunPeriapsis'].values

sunApoapsisTBP = data['thirdBodyAcceleration_sunApoapsis'].values
sunApoapsisSRP = data['srpAcceleration_sunApoapsis'].values

latitude = data['latitude'].values
longitude = data['longitude'].values

x = data['xSurface'].values
y = data['ySurface'].values
z = data['zSurface'].values

# colorsToPlot = U
# # Define the points at the centers of the faces:
# y_coords, x_coords = np.unique(y), np.unique(x)
# y_centers, x_centers = [ arr[:-1] + np.diff(arr)/2 for arr in (y_coords, x_coords)]

# # Convert back to a 2D grid, required for plot_surface:
# print x_coords.size
# print x.shape
# X = x.reshape(-1, x_coords.size)
# Y = y.reshape((y_coords.size, -1))
# Z = z.reshape(X.shape)
# C = colors.reshape(X.shape)
# # the kx, ky define the order of interpolation. Keep it simple, use linear interpolation.
# interp_func = RectBivariateSpline(x_coords, y_coords, C.T, kx=1, ky=1)

# r = ax.plot_surface(X,Y,Z,
#     facecolors=cm.hot(interp_func(x_centers, y_centers).T),
#     rstride=100,  cstride=100)

# xyz = {'x': x, 'y': y, 'z': z}

# # put the data into a pandas DataFrame
# df = pd.DataFrame(xyz, index=range(len(xyz['x'])))

# # re-create the 2D-arrays
# x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
# y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
# x2, y2 = np.meshgrid(x1, y1)
# z2 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')

## Set up the figure for 3D plot
fig = plt.figure( )
ax1 = fig.gca( projection = '3d' )

## Plot the ellipsoidal shape of the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

# latlong = {'longitude':longitude, 'latitude':latitude}
# df = pd.DataFrame(latlong, index=range(len(latlong['longitude'])))
# u = np.linspace(df['longitude'].min(), df['longitude'].max(), len(df['longitude'].unique()))
# v = np.linspace(df['latitude'].min(), df['latitude'].max(), len(df['latitude'].unique()))

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

print U.shape
print ellipsoid_z.shape

# newColor = colors.cnames["slategray"]
# surf = ax1.plot_surface( x2, y2, z2,
#                          rstride=100, cstride=100, antialiased=False )

surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z, rstride=5, cstride=5 )

## format axis and title
ax1.set_xlim( [ -20000, 20000 ] )
ax1.set_ylim( [ -20000, 20000 ] )
ax1.set_zlim( [ -20000, 20000 ] )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)
# azimuthAngle = (Wz * 26358700.0) * 180.0 / np.pi
# print str(azimuthAngle%360)
# ax1.view_init(elev=0, azim=azimuthAngle)

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
