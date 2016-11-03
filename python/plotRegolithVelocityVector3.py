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
import math
from scipy.interpolate import griddata
from numpy import linalg as LA

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

print ""
print "---------------------------------------------------------------------------------"
print "                                 NAOS                                            "
print "                                                                                 "
print "         Copyright (c) 2016, A. Agrawal (abhishek.agrawal@protonmail.com)        "
print "---------------------------------------------------------------------------------"
print ""

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


# Start timer.
start_time = time.time( )

## Set up the figure
fig = plt.figure()
ax1 = fig.gca( projection = '3d' )

# Read data in csv file. data returned as a panda series.
data = pd.read_csv( '../data/singleRegolithEjectaURESolution.csv' )
x = data[ 'x' ].values
y = data[ 'y' ].values
z = data[ 'z' ].values
vx = data[ 'vx' ].values
vy = data[ 'vy' ].values
vz = data[ 'vz' ].values
t = data[ 't' ].values

## Plot the ellipsoidal shape of the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

# newColor = colors.cnames["slategray"]
# surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
#                          rstride=5, cstride=5, alpha=0.5 )
# surf.set_facecolor( newColor )
# surf.set_linewidth( 0.1 )

# ax1.hold( True )

## draw the initial velocity vector
velocityVector = [ vx[0], vy[0], vz[0] ]
positionVector = [ x[0], y[0], z[0] ]

print positionVector

# position sanity check
# positionCheck = positionVector[ 0 ]**2 / alpha**2 + positionVector[ 1 ]**2 / beta**2 + positionVector[ 2 ]**2 / gamma**2
# print positionCheck

# draw the position vector
radius = np.sqrt( ( positionVector[ 0 ] - 0.0 )**2
                + ( positionVector[ 1 ] - 0.0 )**2
                + ( positionVector[ 2 ] - 0.0 )**2 )
ax1.quiver3D( 0.0, 0.0, 0.0,
              positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              length=radius, lw = 1, pivot='tail', arrow_length_ratio=0.2,
              color=colors.cnames["black"], linestyles='solid' )

# draw the north pole vector
ax1.quiver3D( 0.0, 0.0, 0.0,
              0.0, 0.0, gamma,
              length=gamma+1000.0, lw=1, pivot='tail', arrow_length_ratio=0.2,
              color=colors.cnames["black"], linestyles='solid' )

# get the surface unit normal vector (body fixed frame) and plot it
# this is also the z basis vector for the surface frame in body fixed coordinates
normalVector = [ 2.0 * positionVector[ 0 ] / alpha**2,
                 2.0 * positionVector[ 1 ] / beta**2,
                 2.0 * positionVector[ 2 ] / gamma**2 ]

normalVectorMagnitude = LA.norm( normalVector )

unitNormalVector = normalVector / normalVectorMagnitude

ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              unitNormalVector[ 0 ], unitNormalVector[ 1 ], unitNormalVector[ 2 ],
              length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
              color=colors.cnames["darkorange"], linestyles='solid' )

# unitNormal_start = [ positionVector[0], positionVector[1], positionVector[2] ]
# unitNormal_end = [ positionVector[0] + 3000 * unitNormalVector[0],
#                    positionVector[1] + 3000 * unitNormalVector[1],
#                    positionVector[2] + 3000 * unitNormalVector[2] ]
# unitNormal_vecs = list(zip(unitNormal_start, unitNormal_end))
# unitNormal_arrow = Arrow3D(unitNormal_vecs[0],unitNormal_vecs[1],unitNormal_vecs[2], mutation_scale=20, lw=1, arrowstyle="-|>", color="g")
# ax1.add_artist(unitNormal_arrow)

# draw the tangent plane (wireframe plot)
planeConstant = ( 2 * positionVector[ 0 ] ) / alpha**2                                            \
                    + ( 2 * positionVector[ 1 ] ) / beta**2                                       \
                    + ( 2 * positionVector[ 2 ] ) / gamma**2

plane_x = np.linspace( positionVector[ 0 ]-1000, positionVector[ 0 ] + 1000, 10 )
plane_y = np.linspace( positionVector[ 1 ]-1000, positionVector[ 1 ] + 1000, 10 )
planeX, planeY = np.meshgrid( plane_x, plane_y )
planeZ = ( gamma**2 / 2.0 ) * ( planeConstant - ( 2.0 * planeX / alpha**2 ) - ( 2.0 * planeY / beta**2 ) )
ax1.plot_wireframe( planeX, planeY, planeZ, alpha=0.5 )

# get the y basis vector for the surface fixed frame in body frame coordinates
bodyFramePrincipalZAxisUnitVector = [ 0.0, 0.0, 1.0 ]
surfaceFrameYAxisVector = np.cross( unitNormalVector, bodyFramePrincipalZAxisUnitVector )
surfaceFrameYAxisVectorMagnitude = LA.norm( surfaceFrameYAxisVector )
surfaceFrameYAxisUnitVector = surfaceFrameYAxisVector / surfaceFrameYAxisVectorMagnitude

# ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
#               surfaceFrameYAxisUnitVector[ 0 ], surfaceFrameYAxisUnitVector[ 1 ], surfaceFrameYAxisUnitVector[ 2 ],
#               length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
#               color=colors.cnames["green"], linestyles='solid' )

# similarly, now get the x basis vector
surfaceFrameXAxisVector = np.cross( surfaceFrameYAxisUnitVector, unitNormalVector )
surfaceFrameXAxisVectorMagnitude = LA.norm( surfaceFrameXAxisVector )
surfaceFrameXAxisUnitVector = surfaceFrameXAxisVector / surfaceFrameXAxisVectorMagnitude

# ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
#               surfaceFrameXAxisUnitVector[ 0 ], surfaceFrameXAxisUnitVector[ 1 ], surfaceFrameXAxisUnitVector[ 2 ],
#               length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
#               color=colors.cnames["darkred"], linestyles='solid' )

# get the intermediate RTN frame
unitR = positionVector / LA.norm( positionVector )

unitT = np.cross( unitR, bodyFramePrincipalZAxisUnitVector ) / LA.norm( np.cross( unitR, bodyFramePrincipalZAxisUnitVector ) )

unitN = np.cross( unitT, unitR ) / LA.norm( np.cross( unitT, unitR ) )

# ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
#               unitN[ 0 ], unitN[ 1 ], unitN[ 2 ],
#               length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
#               color=colors.cnames["purple"], linestyles='solid' )

# get the angle between the unitN vector (pointing to north but not in tangential plane) and
# the surface normal vector (unitNormalVector)
angle = np.arccos(np.clip(np.dot(unitNormalVector, unitN), -1.0, 1.0))
print angle * 180.0 / np.pi

# cross product of normal vector and unitN
vector1 = np.cross( unitN, unitNormalVector ) / LA.norm( np.cross( unitN, unitNormalVector ) )
# ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
#               vector1[ 0 ], vector1[ 1 ], vector1[ 2 ],
#               length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
#               color=colors.cnames["purple"], linestyles='solid' )

# cross product of normal and unitT
vector2 = np.cross( unitT, unitNormalVector ) / LA.norm( np.cross( unitT, unitNormalVector ) )
ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              vector2[ 0 ], vector2[ 1 ], vector2[ 2 ],
              length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
              color=colors.cnames["purple"], linestyles='solid' )

# cross product of vector2 and normal
vector3 = np.cross( vector2, unitNormalVector ) / LA.norm( np.cross( vector2, unitNormalVector ) )
ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              vector3[ 0 ], vector3[ 1 ], vector3[ 2 ],
              length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
              color=colors.cnames["green"], linestyles='solid' )

# plot the velocity vector now
# ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
#               velocityVector[ 0 ], velocityVector[ 1 ], velocityVector[ 2 ],
#               length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
#               color=colors.cnames["green"], linestyles='solid' )

# vel_start = [ positionVector[0], positionVector[1], positionVector[2] ]
# vel_end = [ positionVector[0] + 500 * velocityVector[0],
#             positionVector[1] + 500 * velocityVector[1],
#             positionVector[2] + 500 * velocityVector[2] ]
# vel_vecs = list(zip(vel_start, vel_end))
# vel_arrow = Arrow3D(vel_vecs[0],vel_vecs[1],vel_vecs[2], mutation_scale=20, lw=1, arrowstyle="-|>", color="g")
# ax1.add_artist(vel_arrow)

## format axis and title
ax1.set_xlim( [ -20000, 20000 ] )
ax1.set_ylim( [ -20000, 20000 ] )
ax1.set_zlim( [ -20000, 20000 ] )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)

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
