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
# from mayavi import mlab
import mayavi.mlab as mLab
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
alpha = 20000.0
beta = 7000.0
gamma = 7000.0

## mayavi 3D plot test
# u = np.linspace(0, 2 * np.pi, 100)
# v = np.linspace(0, np.pi, 100)
# ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
# ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
# ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

# mLab.figure( 1 )
# mLab.mesh( ellipsoid_x, ellipsoid_y, ellipsoid_z, name='asteroid', color=(0.7, 0.5, 0.1) )
# mLab.view( azimuth = -45, elevation = 45 )

# mLab.figure( 2 )
# mLab.mesh( ellipsoid_x, ellipsoid_y, ellipsoid_z, name='asteroid', color=(0.7, 0.5, 0.1) )
# mLab.view( azimuth = -45, elevation = 55 )

# mLab.show( )
# sys.exit( )

## Set up the figure
fig = plt.figure( )
gs = gridspec.GridSpec( 2, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
plt.suptitle( 'Regolith re-impact locations for different launch velocities' )

## Operations
# Connect to SQLite database.
try:
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                    + "multiple_launch_velocity/"
                                    + "simulation_time_9_months/"
                                    + "longestEdge.db")
        # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
        #                            + "multiple_launch_velocity_with_perturbations/"
        #                            + "simulation_time_9_months/"
        #                            + "3.2Density_1cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

azimuthRange = np.arange( 0.0, 360.0, 5.0 )
azimuthRange = tuple( azimuthRange.tolist( ) )

data1 = pd.read_sql( "SELECT    initial_position_x,                                         \
                                initial_position_y,                                         \
                                initial_position_z,                                         \
                                position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                velocity_x,                                                 \
                                velocity_y,                                                 \
                                velocity_z,                                                 \
                                ROUND( initial_velocity_magnitude ),                        \
                                time,                                                       \
                                crash_flag,                                                 \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( crash_flag = 1 )                                          \
                     AND        ROUND( launch_azimuth ) IN" + str(azimuthRange) + ";",      \
                     database )

# data1 = pd.read_sql( "SELECT    initial_position_x,                                         \
#                                 initial_position_y,                                         \
#                                 initial_position_z,                                         \
#                                 position_x,                                                 \
#                                 position_y,                                                 \
#                                 position_z,                                                 \
#                                 velocity_x,                                                 \
#                                 velocity_y,                                                 \
#                                 velocity_z,                                                 \
#                                 ROUND( initial_velocity_magnitude ),                        \
#                                 time,                                                       \
#                                 crash_flag,                                                 \
#                                 ROUND( launch_azimuth )                                     \
#                      FROM       regolith_trajectory_results                                 \
#                      WHERE      ( crash_flag = 1 )                                          \
#                      AND        ROUND( initial_solar_phase_angle ) = 315.0;",               \
#                      database )

data1.columns = [ 'init_pos_x',                                          \
                  'init_pos_y',                                          \
                  'init_pos_z',                                          \
                  'x',                                                   \
                  'y',                                                   \
                  'z',                                                   \
                  'vx',                                                  \
                  'vy',                                                  \
                  'vz',                                                  \
                  'init_vel_mag',                                        \
                  'time',                                                \
                  'crash_flag',                                          \
                  'launch_azimuth' ]

x                   = data1[ 'x' ]
y                   = data1[ 'y' ]
z                   = data1[ 'z' ]
vx                  = data1[ 'vx' ]
vy                  = data1[ 'vy' ]
vz                  = data1[ 'vz' ]
initial_velocity    = data1[ 'init_vel_mag' ]
t                   = data1[ 'time' ]
crashFlag           = data1[ 'crash_flag' ]
azimuth             = data1[ 'launch_azimuth' ]

xPositionStart      = data1[ 'init_pos_x' ]
yPositionStart      = data1[ 'init_pos_y' ]
zPositionStart      = data1[ 'init_pos_z' ]

if database:
    database.close( )

## calculate the lat long for starting point
startRadialDistance = np.sqrt( xPositionStart**2 + yPositionStart**2 + zPositionStart**2 )
startLongitude = np.arctan2( yPositionStart, xPositionStart ) * 180.0 / np.pi
startLatitude = np.arcsin( zPositionStart / startRadialDistance ) * 180.0 / np.pi

## calculate lat long for end points
endRadialDistance = np.sqrt( x**2 + y**2 + z**2 )
endLongitude = np.arctan2( y, x ) * 180 / np.pi
endLatitude = np.arcsin( z / endRadialDistance ) * 180 / np.pi

## indicate starting point
ax1.scatter( startLongitude[0], startLatitude[0], color='black' )
ax1.text( startLongitude[0], startLatitude[0], 'start', size=20, zorder=1, color='black' )

## Plot end points on lat long map
# ax1.scatter( endLongitude, endLatitude, color='red' )
unique_reimpact_velocities = np.unique( initial_velocity )
gridsize = len( unique_reimpact_velocities )

# print gridsize
# print len(initial_velocity)
# print len(endLatitude)
# print len(endLongitude)
# sys.exit( )

hexBinPlot = ax1.hexbin( endLongitude, endLatitude,                                   \
                         C=initial_velocity,                                          \
                         cmap='jet', gridsize=gridsize )
cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax1 )

## format axis and title
ax1.set_xlabel('longitude [degree]')
ax1.set_ylabel('latitude [degree]')
cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )
ax1.set_xlim( -180.0, 180.0 )
ax1.set_yticks( np.arange( -90.0, 90.0, 10.0 ) )
ax1.set_title( 'Qualitative re-impact map' )
ax1.grid( True )

## plot the actual data points, scatter points
# ax2.scatter( endLongitude, endLatitude, s=5, c=initial_velocity, edgecolors='face', cmap='jet' )
# cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax2 )

unique_reimpact_velocities = np.unique( initial_velocity ).tolist( )
print "Final re-impact velocity = " + str( unique_reimpact_velocities[ -1 ] )
pltColors = plt.cm.Vega20( np.linspace( 0, 1, len( unique_reimpact_velocities ) ) )

for index in range( 0, len( unique_reimpact_velocities ) ):
    current_velocity = unique_reimpact_velocities[ index ]
    current_velocity_indices = np.where( initial_velocity == current_velocity )
    current_velocity_indices = current_velocity_indices[ 0 ]
    plotLatitude = endLatitude[ current_velocity_indices ]
    plotLongitude = endLongitude[ current_velocity_indices ]
    ax2.scatter( plotLongitude, plotLatitude, s=5, c=pltColors[ index ],                    \
                 edgecolors='face',                                                         \
                 label='$V_{launch}$ = ' + str( current_velocity ) + ' [m/s]' )

ax2.legend( markerscale=7 ).draggable( )
ax2.set_xlabel('longitude [degree]')
ax2.set_ylabel('latitude [degree]')
# cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )
ax2.set_xlim( -180.0, 180.0 )
ax2.set_yticks( np.arange( -90.0, 90.0, 10.0 ) )
ax2.set_title( 'Actual re-impact locations' )
ax2.grid( True )

## 2D reimpact locations for specific velocities
fig = plt.figure( figsize=(10, 8) )
gs = gridspec.GridSpec( 3, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
plt.suptitle( "Regolith re-impact behavior with varying launch azimuth" )

scatterVelocity = [ 5.0, 9.0, 13.0 ]
plotHandles = [ ax1, ax2, ax3 ]
for index in range( 0, len( scatterVelocity ) ):
    currentVelocity = scatterVelocity[ index ]
    currentPlotHandle = plotHandles[ index ]
    data_indices = np.where( initial_velocity == currentVelocity )
    data_indices = data_indices[ 0 ]

    plotAzimuth = azimuth[ data_indices ]
    plotLatitude = endLatitude[ data_indices ]
    plotLongitude = endLongitude[ data_indices ]

    singleVelocityScatter = currentPlotHandle.scatter( plotLongitude, plotLatitude,
                                         c=plotAzimuth,
                                         cmap='jet' )
    cbar = plt.colorbar( singleVelocityScatter, cmap='jet', ax=currentPlotHandle )

    currentPlotHandle.set_title( 'Launch velocity = ' + str( currentVelocity ) + ' [m/s]' )
    currentPlotHandle.set_xlabel('longitude [degree]')
    currentPlotHandle.set_ylabel('latitude [degree]')
    cbar.ax.set_ylabel( 'Launch azimuth [degree]' )
    currentPlotHandle.set_xlim( -180.0, 180.0 )
    currentPlotHandle.set_yticks( np.arange( -90.0, 90.0, 20.0 ) )
    currentPlotHandle.grid( True )

plt.show( )
sys.exit( )

## 3D crash point depiction
# Set up the figure
fig = plt.figure( figsize=(12, 10) )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ], projection='3d' )
ax2 = plt.subplot( gs[ 1 ], projection='3d' )
ax3 = plt.subplot( gs[ 2 ], projection='3d' )
ax4 = plt.subplot( gs[ 3 ], projection='3d' )
plt.suptitle('Re-impact locations on ellipsoidal asteroid')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))
newColor = colors.cnames["lightgray"]
plotVelMag = [5.0, 9.0, 13.0]

plotHandles = [ ax1, ax2, ax3, ax4 ]
for index in range( 0, len( plotHandles ) ):
    currentPlotHandle = plotHandles[ index ]
    surf1 = currentPlotHandle.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                                            zorder=10,
                                            rstride=5, cstride=5, alpha=1.0,
                                            color=newColor )
    surf1.set_linewidth( 1.0 )
    currentPlotHandle.scatter( xPositionStart[ 0 ], yPositionStart[ 0 ], zPositionStart[ 0 ],
                 zorder=100,
                 s=5,
                 c='black',
                 label='Launch site' )

    for velIndex in range( 0, len( plotVelMag ) ):
        currentVelocity = plotVelMag[ velIndex ]
        data_indices = np.where( initial_velocity == currentVelocity )
        data_indices = data_indices[ 0 ]
        crash_position_x = x[data_indices]
        crash_position_y = y[data_indices]
        crash_position_z = z[data_indices]

        currentPlotHandle.scatter( crash_position_x, crash_position_y, crash_position_z,
                     zorder=10,
                     s=5,
                     c=pltColors[ int(currentVelocity) - 1 ],
                     label='$V_{launch}$ = ' + str( currentVelocity ) + ' [m/s]' )

        currentPlotHandle.set_xlabel('x [m]')
        currentPlotHandle.set_ylabel('y [m]')
        currentPlotHandle.set_zlabel('z [m]')
        # currentPlotHandle.set_title('Re-impact locations for regolith')
        currentPlotHandle.set_xlim( -alpha, alpha )
        currentPlotHandle.set_ylim( -alpha, alpha )
        currentPlotHandle.set_zlim( -alpha, alpha )
        currentPlotHandle.legend( markerscale=7 ).draggable( )

## MAYAVI plots
# plotVelocity = 5.0
# data_indices = np.where( initial_velocity == plotVelocity )
# data_indices = data_indices[ 0 ]
# crash_x = x[ data_indices ]
# crash_y = y[ data_indices ]
# crash_z = z[ data_indices ]
# mLab.mesh( ellipsoid_x, ellipsoid_y, ellipsoid_z )
# mLab.points3d( crash_x, crash_y, crash_z )

plotVelMag = [5.0, 9.0, 13.0]
viewAzimuth = [ 45, 135, 225, 315 ]
regolithColor = [ (1.0, 0, 0), (0, 1.0, 0), (0, 0, 1.0) ]
for velIndex in range( 0, len( plotVelMag ) ):
    currentRegolithColor = regolithColor[ velIndex ]
    currentVelocity = plotVelMag[ velIndex ]
    data_indices = np.where( initial_velocity == currentVelocity )
    data_indices = data_indices[ 0 ]
    crash_position_x = x[data_indices]
    crash_position_y = y[data_indices]
    crash_position_z = z[data_indices]

    for plotIndex in range( 0, len( viewAzimuth ) ):
        mLab.figure( bgcolor = (0, 0, 0) )
        mLab.mesh( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                   color=(140.0/255.0, 140.0/255.0, 140.0/255.0) )
        mLab.points3d( xPositionStart[ 0 ], yPositionStart[ 0 ], zPositionStart[ 0 ],
                       scale_factor=500,
                       color=(0, 0, 0) )
        mLab.points3d( crash_position_x, crash_position_y, crash_position_z,
                       scale_factor=500,
                       color=currentRegolithColor )
        mLab.view( azimuth = viewAzimuth[ plotIndex ], elevation = 75 )

    mLab.show( )

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
