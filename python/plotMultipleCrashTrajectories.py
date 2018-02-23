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

def plotEllipse( semiMajor, semiMinor, plotHandle ):
        angleRange = np.linspace(0, 2 * np.pi, 360)
        r = semiMajor * semiMinor / np.sqrt( ( semiMinor * np.cos( angleRange ) )**2
                                           + ( semiMajor * np.sin( angleRange ) )**2 )
        x = r * np.cos( angleRange )
        y = r * np.sin( angleRange )
        plotHandle.plot( x, y, label='Asteroid at launch', linestyle='dotted', color='black', lw=0.5 )
        plotHandle.grid( )

# Start timer.
start_time = time.time( )

## Operations
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
rotationRate = 0.00033118202125129593

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

## Ellipse patches example
# fig = plt.figure( )
# gs = gridspec.GridSpec( 1, 1 )
# ax1 = plt.subplot( gs[ 0 ] )
# ellipse = mpatches.Ellipse( ( 0.0, 0.0 ), 40.0e3, 14.0e3, 135.0,
#                             facecolor='none',
#                             fill=False,
#                             edgecolor='black',
#                             lw=2.,
#                             ls='solid',
#                             label='Asteroid at re-impact' )
# ax1.add_patch( ellipse )
# ax1.plot( [], [] )
# ax1.set_aspect( 1 )
# plt.show( )
# sys.exit( )

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
        #                            + "3.2Density_1cmSize/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

# azimuthRange = np.arange(1.0, 91.0, 1.0)
# azimuthRange = tuple( azimuthRange.tolist( ) )
azimuthRange = (270.0)

print "Extracting data now...\n"

data = pd.read_sql( "SELECT     trajectory_id,                                              \
                                ROUND( launch_azimuth ),                                    \
                                time                                                        \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ROUND( crash_flag ) = 1                                     \
                     AND        ROUND( initial_velocity_magnitude ) = 5.0                   \
                     AND        ROUND( launch_azimuth ) IN " + str( azimuthRange ) + ";",   \
                     database )

data.columns = [ 'trajectory_id', 'launch_azimuth', 'time' ]

trajectory_id = data[ 'trajectory_id' ]
finalTime     = data[ 'time' ]
launchAzim    = data[ 'launch_azimuth' ]
trajectory_id = tuple( trajectory_id.tolist( ) )

data1 = pd.read_sql( "SELECT    initial_position_x,                                         \
                                initial_position_y,                                         \
                                initial_position_z,                                         \
                                inertial_position_x,                                        \
                                inertial_position_y,                                        \
                                inertial_position_z,                                        \
                                position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                ROUND( initial_velocity_magnitude ),                        \
                                time,                                                       \
                                crash_flag,                                                 \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      trajectory_id IN " + str( trajectory_id ) + ";",            \
                     database )

data1.columns = [ 'init_pos_x',                                          \
                  'init_pos_y',                                          \
                  'init_pos_z',                                          \
                  'x',                                                   \
                  'y',                                                   \
                  'z',                                                   \
                  'bodyFrame_x',                                         \
                  'bodyFrame_y',                                         \
                  'bodyFrame_z',                                         \
                  'init_vel_mag',                                        \
                  'time',                                                \
                  'crash_flag',                                          \
                  'launch_azimuth' ]

x                   = data1[ 'x' ]
y                   = data1[ 'y' ]
z                   = data1[ 'z' ]
bodyFrame_x         = data1[ 'bodyFrame_x' ]
bodyFrame_y         = data1[ 'bodyFrame_y' ]
bodyFrame_z         = data1[ 'bodyFrame_z' ]
initial_velocity    = data1[ 'init_vel_mag' ]
t                   = data1[ 'time' ]
crashFlag           = data1[ 'crash_flag' ]
azimuth             = data1[ 'launch_azimuth' ]

xPositionStart      = data1[ 'init_pos_x' ]
yPositionStart      = data1[ 'init_pos_y' ]
zPositionStart      = data1[ 'init_pos_z' ]

if database:
    database.close( )

# Extraction time
end_time = time.time( )
print "Data extraction time: " + str("{:,g}".format(end_time - start_time)) + "s\n"

print "Plotting data now...\n"

viewAngles = [ 45 ] # for the mayavi 3D view

## Plot all crash trajectories for the given azimuth angle range and launch velocity
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

uniqueAzimuths = np.unique( azimuth )
for azimuthIndex in range( 0, len( uniqueAzimuths ) ):
    currentAzimuth = uniqueAzimuths[ azimuthIndex ]
    data_indices = np.where( azimuth == currentAzimuth )
    data_indices = data_indices[ 0 ]
    trajectory_x = x[data_indices]
    trajectory_y = y[data_indices]
    trajectory_z = z[data_indices]

    ax1.plot( trajectory_x, trajectory_y,
              color='red',
              linewidth=0.5 )

## Plot the trajectory which took the maximum amount of time to crash within the given azimuth range
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

maxTimeValue = max(t)

for azimuthIndex in range( 0, len( uniqueAzimuths ) ):
    currentAzimuth = uniqueAzimuths[ azimuthIndex ]
    data_indices = np.where( azimuth == currentAzimuth )
    data_indices = data_indices[ 0 ]
    trajectory_x = x[data_indices]
    trajectory_y = y[data_indices]
    trajectory_z = z[data_indices]

    currentTimeValues = t[data_indices]
    currentTimeValues = currentTimeValues.tolist( )
    currentLastTimeValue = currentTimeValues[ -1 ]

    if currentLastTimeValue == maxTimeValue:
        plotEllipse( alpha, beta, ax1 )
        plotAzimuth = azimuth[data_indices]

        # calculate lat long for end points
        dataIndex = np.where( t == currentLastTimeValue )
        dataIndex = dataIndex[ 0 ]
        body_x = bodyFrame_x[ dataIndex ]
        body_y = bodyFrame_y[ dataIndex ]
        body_z = bodyFrame_z[ dataIndex ]
        endRadialDistance = np.sqrt( body_x**2 + body_y**2 + body_z**2 )
        endLongitude = np.arctan2( body_y, body_x ) * 180 / np.pi
        endLatitude = np.arcsin( body_z / endRadialDistance ) * 180 / np.pi

        print "Longitude = " + str(endLongitude)
        print "Latitude = " + str(endLatitude)

        ax1.plot( trajectory_x, trajectory_y,
                  color='red',
                  linewidth=1.0,
                  label='Launzh azimuth = ' + str(currentAzimuth) + ' [deg]' )

        acwRotationAngle = rotationRate * currentLastTimeValue * 180.0 / np.pi
        ellipse = mpatches.Ellipse( ( 0.0, 0.0 ), 40.0e3, 14.0e3, acwRotationAngle, alpha=1.0,
                                    fill=False,
                                    edgecolor='black',
                                    label='Asteroid at re-impact' )
        ax1.add_patch( ellipse )
        ax1.plot( [], [] )

        ax1.set_aspect( 1 )
        ax1.grid(True)
        ax1.legend( ).draggable( )
        ax1.set_xlabel('x [m]')
        ax1.set_ylabel('y [m]')
        ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)

        # mLab.figure( bgcolor = (0, 0, 0) )
        # mLab.mesh( ellipsoid_x, ellipsoid_y, ellipsoid_z,
        #            color=(140.0/255.0, 140.0/255.0, 140.0/255.0) )
        # mLab.plot3d( trajectory_x, trajectory_y, trajectory_z,
        #              color=(1.0, 0.0, 0.0),
        #              line_width=1.5,
        #              representation='points' )

## Plot the trajectory which took the least amount of time to crash within the given azimuth range
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

minTimeValue = min(finalTime)
minTimeValue = minTimeValue.tolist( )
print '\n'
print minTimeValue
print '\n'

dataIndex = np.where( t == minTimeValue )
dataIndex = dataIndex[ 0 ]
minTimeAzimuth = azimuth[ dataIndex ]
minTimeAzimuth = minTimeAzimuth.tolist( )

body_x = bodyFrame_x[ dataIndex ]
body_y = bodyFrame_y[ dataIndex ]
body_z = bodyFrame_z[ dataIndex ]
endRadialDistance = np.sqrt( body_x**2 + body_y**2 + body_z**2 )
endLongitude = np.arctan2( body_y, body_x ) * 180 / np.pi
endLatitude = np.arcsin( body_z / endRadialDistance ) * 180 / np.pi

print "Longitude = " + str(endLongitude)
print "Latitude = " + str(endLatitude)

data_indices = np.where( azimuth == minTimeAzimuth[ 0 ] )
data_indices = data_indices[ 0 ]
trajectory_x = x[data_indices]
trajectory_y = y[data_indices]
trajectory_z = z[data_indices]

plotEllipse( alpha, beta, ax1 )

ax1.plot( trajectory_x, trajectory_y,
          color='red',
          linewidth=1.0,
          label='Launzh azimuth = ' + str(minTimeAzimuth) + ' [deg]' )

acwRotationAngle = rotationRate * minTimeValue * 180.0 / np.pi
ellipse2 = mpatches.Ellipse( ( 0.0, 0.0 ), 40.0e3, 14.0e3, acwRotationAngle, alpha=1.0,
                            fill=False,
                            edgecolor='black',
                            label='Asteroid at re-impact' )
ax1.add_patch( ellipse2 )
ax1.plot( [], [] )

ax1.set_aspect( 1 )
ax1.grid(True)
ax1.legend( ).draggable( )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)

## plot all total time values against launch azimuth for the given azimuth angle range
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

ax1.scatter( launchAzim, finalTime, s=7, c='black' )
ax1.grid(True)
ax1.set_xlabel('Launch azimuth [deg]')
ax1.set_ylabel('Time to re-impact [s]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## Show the plot
plt.show( )
# mLab.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
