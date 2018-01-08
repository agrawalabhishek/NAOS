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

def plotEllipse( semiMajor, semiMinor, plotHandle ):
        angleRange = np.linspace(0, 2 * np.pi, 360)
        r = semiMajor * semiMinor / np.sqrt( ( semiMinor * np.cos( angleRange ) )**2
                                           + ( semiMajor * np.sin( angleRange ) )**2 )
        x = r * np.cos( angleRange )
        y = r * np.sin( angleRange )
        plotHandle.plot( x, y, label='Ellipsoid 2D projection' )
        plotHandle.grid( )

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

# Connect to SQLite database.
try:
    # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
    #                            + "multiple_launch_velocity_with_perturbations/"
    #                            + "simulation_time_9_months/"
    #                            + "test.db")
    database = sqlite3.connect("../data/VandV/"
                                + "new/"
                                + "launch_velocity/"
                                + "leadingEdge2.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

lowerTime = 0.0
upperTime = 270.0 * 24.0 * 60.0 * 60.0
solarPhase = 315.0

data = pd.read_sql( "SELECT     position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                velocity_x,                                                     \
                                velocity_y,                                                     \
                                velocity_z,                                                     \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( initial_solar_phase_angle ) = " + str(solarPhase) + "    \
                     AND        time >= " + str(lowerTime)
                                + " AND time <= " + str(upperTime) + "                          \
                     AND        ROUND( launch_azimuth ) = 135.0                                 \
                     AND        ROUND( launch_declination ) = 60.0                              \
                     AND        ROUND( initial_velocity_magnitude ) = 6.0;",                    \
                     database )

data.columns = [ 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'velocity_x',                                          \
                 'velocity_y',                                          \
                 'velocity_z',                                          \
                 'time' ]

if database:
    database.close( )

x                   = data[ 'x' ]
y                   = data[ 'y' ]
z                   = data[ 'z' ]
vx                  = data[ 'velocity_x' ]
vy                  = data[ 'velocity_y' ]
vz                  = data[ 'velocity_z' ]
t                   = data[ 'time' ]

print "Processing data now...\n"

## Plot the ellipsoidal shape of the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

newColor = colors.cnames["slategray"]
surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         color=newColor,
                         rstride=5, cstride=5, alpha=0.2 )
# surf.set_facecolor( newColor )
surf.set_linewidth( 3.0 )

## draw the initial velocity vector
velocityVector = [ vx[0], vy[0], vz[0] ]
positionVector = [ x[0], y[0], z[0] ]

# # draw the position vector
# radius = np.sqrt( ( positionVector[ 0 ] - 0.0 )**2
#                 + ( positionVector[ 1 ] - 0.0 )**2
#                 + ( positionVector[ 2 ] - 0.0 )**2 )
# ax1.quiver3D( 0.0, 0.0, 0.0,
#               positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
#               length=radius, lw = 1, pivot='tail', arrow_length_ratio=0.2,
#               color=colors.cnames["black"], linestyles='solid' )

# # draw the north pole vector
# ax1.quiver3D( 0.0, 0.0, 0.0,
#               0.0, 0.0, gamma,
#               length=500, lw=1, pivot='tail', arrow_length_ratio=0.01,
#               color=colors.cnames["black"], linestyles='solid' )
ax1.quiver3D( 0.0, 0.0, 0.0,
              0.0, 0.0, gamma,
              lw=1, pivot='tail',
              color=colors.cnames["black"], linestyles='solid', label='North Pole' )

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
              color=colors.cnames["darkorange"], linestyles='solid', label='Surface Normal / Z-axis local' )

# get the intermediate RTN frame
unitR = positionVector / LA.norm( positionVector )

bodyFramePrincipalZAxisUnitVector = [ 0.0, 0.0, 1.0 ]
unitT = np.cross( unitR, bodyFramePrincipalZAxisUnitVector ) / LA.norm( np.cross( unitR, bodyFramePrincipalZAxisUnitVector ) )

unitN = np.cross( unitT, unitR ) / LA.norm( np.cross( unitT, unitR ) )

# check for if the position vector is pointing to the poles
unitPositionVector = positionVector / LA.norm( positionVector )
positionDotPrincipalZ = np.dot( unitPositionVector, [ 0.0, 0.0, 1.0 ] )
positionDotNegativePrincipalZ = np.dot( unitPositionVector, [ 0.0, 0.0, -1.0 ] )

# cross product of normal and unitT
unitX = np.cross( unitT, unitNormalVector ) / LA.norm( np.cross( unitT, unitNormalVector ) )

# cross product of unitX and normal
unitY = np.cross( unitNormalVector, unitX ) / LA.norm( np.cross( unitNormalVector, unitX ) )

if positionDotPrincipalZ == 1.0:
    unitX = [ 1.0, 0.0, 0.0 ]
    unitY = [ 0.0, 1.0, 0.0 ]
elif positionDotNegativePrincipalZ == 1.0:
    unitX = [ -1.0, 0.0, 0.0 ]
    unitY = [ 0.0, 1.0, 0.0 ]

#plot the x and y basis vectors
ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              unitX[ 0 ], unitX[ 1 ], unitX[ 2 ],
              length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
              color=colors.cnames["magenta"], linestyles='solid', label='Local North direction / X-axis local' )

ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              unitY[ 0 ], unitY[ 1 ], unitY[ 2 ],
              length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
              color=colors.cnames["blue"], linestyles='solid', label='Y-axis local' )

# plot the velocity vector now
ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              velocityVector[ 0 ], velocityVector[ 1 ], velocityVector[ 2 ],
              length=500,
              lw=1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["green"],
              linestyles='solid',
              label='Velocity' )

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
ax1.legend( ).draggable( )

## 2D plots for the RTN frame
fig = plt.figure( figsize=(7,10) )
plt.suptitle( "Local surface frame at longest edge of the asteroid" )
gs  = gridspec.GridSpec( 3, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

# xz plot
plotEllipse( alpha, gamma, ax1 )

ax1.quiver( positionVector[ 0 ], positionVector[ 2 ],
            unitX[ 0 ], unitX[ 2 ],
            # units='xy',
            angles='xy',
            width=0.005,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color=colors.cnames["magenta"],
            label='Local North direction / X-axis local' )

ax1.quiver( positionVector[ 0 ], positionVector[ 2 ],
            unitNormalVector[ 0 ], unitNormalVector[ 2 ],
            # units='xy',
            angles='xy',
            width=0.005,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color=colors.cnames["darkorange"],
            label='Surface Normal / Z-axis local' )

ax1.quiver( positionVector[ 0 ], positionVector[ 2 ],
            unitY[ 0 ], unitY[ 2 ],
            # units='xy',
            angles='xy',
            width=0.005,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color=colors.cnames["blue"],
            label='Y-axis local' )

ax1.grid(True)
# ax1.legend( ).draggable( )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('z [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)

# xy plot
plotEllipse( alpha, beta, ax2 )

ax2.quiver( positionVector[ 0 ], positionVector[ 1 ],
            unitX[ 0 ], unitX[ 1 ],
            # units='xy',
            angles='xy',
            width=0.005,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color=colors.cnames["magenta"],
            label='Local North direction / X-axis local' )

ax2.quiver( positionVector[ 0 ], positionVector[ 1 ],
            unitNormalVector[ 0 ], unitNormalVector[ 1 ],
            # units='xy',
            angles='xy',
            width=0.005,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color=colors.cnames["darkorange"],
            label='Surface Normal / Z-axis local' )

ax2.quiver( positionVector[ 0 ], positionVector[ 1 ],
            unitY[ 0 ], unitY[ 1 ],
            # units='xy',
            angles='xy',
            width=0.005,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color=colors.cnames["blue"],
            label='Y-axis local' )

ax2.grid(True)
ax2.legend( ).draggable( )
ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)

# yz plot
plotEllipse( beta, gamma, ax3 )

ax3.quiver( positionVector[ 1 ], positionVector[ 2 ],
            unitX[ 1 ], unitX[ 2 ],
            # units='xy',
            angles='xy',
            width=0.005,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color=colors.cnames["magenta"],
            label='Local North direction / X-axis local' )

ax3.quiver( positionVector[ 1 ], positionVector[ 2 ],
            unitNormalVector[ 1 ], unitNormalVector[ 2 ],
            # units='xy',
            angles='xy',
            width=0.005,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color=colors.cnames["darkorange"],
            label='Surface Normal / Z-axis local' )

ax3.quiver( positionVector[ 1 ], positionVector[ 2 ],
            unitY[ 1 ], unitY[ 2 ],
            # units='xy',
            angles='xy',
            width=0.005,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color=colors.cnames["blue"],
            label='Y-axis local' )

ax3.grid(True)
# ax3.legend( ).draggable( )
ax3.set_xlabel('y [m]')
ax3.set_ylabel('z [m]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)

# check for normal vector and surface frame x-axis
# works only if its on the longest edge coz then its aligned with the position vector
normalVectorCheck = np.cross( unitNormalVector, unitR )
print normalVectorCheck

xSurfaceFrameCheck = np.dot( unitNormalVector, unitX )
print xSurfaceFrameCheck

# check if normal and local north direction are perpendicular - works for any point on ellipsoid
normalToLocalNorthCheck = np.dot( unitX, unitNormalVector )
print normalToLocalNorthCheck

# XY projection of local north direction vector (surface x axis vector) and the position vector
# will be aligned if the x axis vector points in the north
northDirectionCheckGeneral = [ unitX[ 0 ]/unitR[ 0 ], unitX[ 1 ]/unitR[ 1 ] ]
print northDirectionCheckGeneral

print unitX
print unitR
print unitNormalVector

## Velocity angle check (Very Important)
fig = plt.figure()
ax1 = fig.gca( projection = '3d' )

print "\n\n"

unitVelocityVector = velocityVector / LA.norm( velocityVector )
declinationAngle = np.arccos(np.clip(np.dot(unitVelocityVector, unitNormalVector), -1.0, 1.0)) * 180 / np.pi

# projection of velocity vector ONTO a PLANE (XY plane of the surface frame)
velocityVector_normalVectorProjection = np.multiply(((np.dot(velocityVector, normalVector)) / LA.norm(normalVector)**2), normalVector)
velocityVector_xyProjection = velocityVector - velocityVector_normalVectorProjection
unitVelocityVector_xyProjection = velocityVector_xyProjection / LA.norm(velocityVector_xyProjection)

azimuthAngle = np.arccos(np.clip(np.dot(unitVelocityVector_xyProjection, unitX), -1.0, 1.0)) * 180 / np.pi

print "Velocity Magnitude = " + str(LA.norm(velocityVector))
print "velocity vector" + str(velocityVector)
print "velocity xy projection = " + str(velocityVector_xyProjection)
print "x basis = " + str(unitX)
print "unit normal vector = " + str(unitNormalVector)
print "Declination = " + str( declinationAngle ) + " [deg]"
print "Azimuth = " + str( azimuthAngle ) + " [deg]"

# plot ellipsoid
# surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
#                          color=newColor,
#                          rstride=5, cstride=5, alpha=0.2 )
# surf.set_linewidth( 3.0 )

# plot the x, y, z basis vectors of the surface frame
ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              unitX[ 0 ], unitX[ 1 ], unitX[ 2 ],
              length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
              color=colors.cnames["magenta"], linestyles='solid', label='Local North direction / X-axis local' )

ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              unitY[ 0 ], unitY[ 1 ], unitY[ 2 ],
              length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
              color=colors.cnames["blue"], linestyles='solid', label='Y-axis local' )

ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              unitNormalVector[ 0 ], unitNormalVector[ 1 ], unitNormalVector[ 2 ],
              length=3000, lw=1, pivot='tail', arrow_length_ratio=0.2,
              color=colors.cnames["darkorange"], linestyles='solid', label='Surface Normal / Z-axis local' )

# plot the velocity vector now
ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              velocityVector[ 0 ], velocityVector[ 1 ], velocityVector[ 2 ],
              length=500,
              lw=1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["green"],
              linestyles='solid',
              label='Velocity' )

ax1.quiver3D( positionVector[ 0 ], positionVector[ 1 ], positionVector[ 2 ],
              velocityVector_xyProjection[ 0 ],
              velocityVector_xyProjection[ 1 ],
              velocityVector_xyProjection[ 2 ],
              length=500,
              lw=1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["black"],
              linestyles='solid',
              label='Velocity vector projection' )

# plot the rotating frame
ax1.quiver3D( 0.0, 0.0, 0.0,
              alpha,
              0.0,
              0.0,
              # length=50,
              lw=1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["red"],
              linestyles='dotted',
              label='x-axis ARF' )

ax1.quiver3D( 0.0, 0.0, 0.0,
              0.0,
              beta,
              0.0,
              # length=50,
              lw=1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["brown"],
              linestyles='dotted',
              label='y-axis ARF' )

ax1.quiver3D( 0.0, 0.0, 0.0,
              0.0,
              0.0,
              gamma,
              # length=50,
              lw=1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["blue"],
              linestyles='dotted',
              label='z-axis ARF' )

# draw position vector to launch point
ax1.quiver3D( 0.0, 0.0, 0.0,
              positionVector[ 0 ],
              positionVector[ 1 ],
              positionVector[ 2 ],
              # length=50,
              lw=1,
              pivot='tail',
              arrow_length_ratio=0.2,
              color=colors.cnames["black"],
              linestyles='dotted',
              label='Launch location' )

ax1.set_xlim( [ -20000, 20000 ] )
ax1.set_ylim( [ -20000, 20000 ] )
ax1.set_zlim( [ -20000, 20000 ] )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), offset=False)
ax1.legend( ).draggable( )

## Show the plot
# plt.show( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
