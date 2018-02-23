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
muAsteroid = 876514

## Operations
# Connect to SQLite database.
try:
        # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
        #                             + "multiple_launch_velocity/"
        #                             + "simulation_time_9_months/"
        #                             + "longestEdge.db")
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "3.2Density_1cmRadius/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

azimuthRange = np.arange( 0.0, 360.0, 5.0 )
azimuthRange = tuple( azimuthRange.tolist( ) )

# data1 = pd.read_sql( "SELECT    initial_position_x,                                         \
#                                 initial_position_y,                                         \
#                                 initial_position_z,                                         \
#                                 position_x,                                                 \
#                                 position_y,                                                 \
#                                 position_z,                                                 \
#                                 velocity_x,                                                 \
#                                 velocity_y,                                                 \
#                                 velocity_z,                                                 \
#                                 initial_inertial_velocity_x,                                \
#                                 initial_inertial_velocity_y,                                \
#                                 initial_inertial_velocity_z,                                \
#                                 directional_inertial_escape_speed,                          \
#                                 sma,                                                        \
#                                 ROUND( initial_velocity_magnitude ),                        \
#                                 time,                                                       \
#                                 ROUND( launch_azimuth )                                     \
#                      FROM       regolith_trajectory_results                                 \
#                      WHERE      ( escape_flag = 1 )                                         \
#                      AND        ROUND( launch_azimuth ) IN" + str(azimuthRange) + ";",      \
#                      database )

data1 = pd.read_sql( "SELECT    initial_position_x,                                         \
                                initial_position_y,                                         \
                                initial_position_z,                                         \
                                position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                velocity_x,                                                 \
                                velocity_y,                                                 \
                                velocity_z,                                                 \
                                initial_inertial_velocity_x,                                \
                                initial_inertial_velocity_y,                                \
                                initial_inertial_velocity_z,                                \
                                directional_inertial_escape_speed,                          \
                                sma,                                                        \
                                ROUND( initial_velocity_magnitude ),                        \
                                time,                                                       \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 )                                         \
                     AND        ROUND( initial_solar_phase_angle ) = 45;",                  \
                     database )

data1.columns = [ 'init_pos_x',                                          \
                  'init_pos_y',                                          \
                  'init_pos_z',                                          \
                  'x',                                                   \
                  'y',                                                   \
                  'z',                                                   \
                  'vx',                                                  \
                  'vy',                                                  \
                  'vz',                                                  \
                  'initial_inertial_velocity_x',                         \
                  'initial_inertial_velocity_y',                         \
                  'initial_inertial_velocity_z',                         \
                  'directional_inertial_escape_speed',                   \
                  'sma',                                                 \
                  'init_vel_mag',                                        \
                  'time',                                                \
                  'launch_azimuth' ]

x                   = data1[ 'x' ]
y                   = data1[ 'y' ]
z                   = data1[ 'z' ]
vx                  = data1[ 'vx' ]
vy                  = data1[ 'vy' ]
vz                  = data1[ 'vz' ]
inertial_vx         = data1[ 'initial_inertial_velocity_x' ]
inertial_vy         = data1[ 'initial_inertial_velocity_y' ]
inertial_vz         = data1[ 'initial_inertial_velocity_z' ]
sma                 = data1[ 'sma' ]
initial_velocity    = data1[ 'init_vel_mag' ]
t                   = data1[ 'time' ]
azimuth             = data1[ 'launch_azimuth' ]

xPositionStart      = data1[ 'init_pos_x' ]
yPositionStart      = data1[ 'init_pos_y' ]
zPositionStart      = data1[ 'init_pos_z' ]

inertial_guaranteedEscape = data1[ 'directional_inertial_escape_speed' ]

if database:
    database.close( )

fig = plt.figure( figsize=(10, 8) )
gs = gridspec.GridSpec( 3, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
plt.suptitle( "Regolith escape behavior with varying launch azimuth" )

data_indices = np.where( initial_velocity == 5.0 )
data_indices = data_indices[ 0 ]
hev = np.sqrt( -muAsteroid / sma[ data_indices ] )
ax1.scatter( azimuth[data_indices], hev )
print str( max(hev/5.0) * 100.0 )

data_indices = np.where( initial_velocity == 9.0 )
data_indices = data_indices[ 0 ]
hev = np.sqrt( -muAsteroid / sma[ data_indices ] )
ax2.scatter( azimuth[data_indices], hev )
print str( max(hev/9.0) * 100.0 )

data_indices = np.where( initial_velocity == 13.0 )
data_indices = data_indices[ 0 ]
hev = np.sqrt( -muAsteroid / sma[ data_indices ] )
ax3.scatter( azimuth[data_indices], hev )
print str( max(hev/13.0) * 100.0 )

fig = plt.figure( figsize=(10, 8) )
gs = gridspec.GridSpec( 3, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

data_indices = np.where( initial_velocity == 5.0 )
data_indices = data_indices[ 0 ]
hev = np.sqrt( -muAsteroid / sma[ data_indices ] )
ax1.hist( hev, 50, facecolor='green' )

data_indices = np.where( initial_velocity == 9.0 )
data_indices = data_indices[ 0 ]
hev = np.sqrt( -muAsteroid / sma[ data_indices ] )
ax2.hist( hev, 50, facecolor='green' )

data_indices = np.where( initial_velocity == 13.0 )
data_indices = data_indices[ 0 ]
hev = np.sqrt( -muAsteroid / sma[ data_indices ] )
ax3.hist( hev, 50, facecolor='green' )

## Combined plots
fig = plt.figure( )
plt.suptitle("Escape characteristic behavior, Launch velocity = 5.0 [m/s]")
data_indices = np.where( initial_velocity == 5.0 )
data_indices = data_indices[ 0 ]
hev = np.sqrt( -muAsteroid / sma[ data_indices ] )

ax1 = plt.subplot2grid((8,15), (0, 0), rowspan=8, colspan=10)
ax1.scatter( azimuth[data_indices], hev )
ax1.grid(True)
ax1.set_xlabel('Launch azimuth [deg]')
ax1.set_ylabel('Hyperbolic excess velocity [m/s]')

ax2 = plt.subplot2grid((8,15), (0, 10), rowspan=8, colspan=5, sharey=ax1)
ax2.hist( hev, 40,
          facecolor=colors.cnames['magenta'], edgecolor='black',
          orientation='horizontal', alpha=1.0 )
plt.setp( ax2.get_yticklabels( ), visible=False )
ax2.grid(True)
ax2.set_xlabel('Regolith count')

fig = plt.figure( )
plt.suptitle("Escape characteristic behavior, Launch velocity = 9.0 [m/s]")
data_indices = np.where( initial_velocity == 9.0 )
data_indices = data_indices[ 0 ]
hev = np.sqrt( -muAsteroid / sma[ data_indices ] )

ax1 = plt.subplot2grid((8,15), (0, 0), rowspan=8, colspan=10)
ax1.scatter( azimuth[data_indices], hev )
ax1.grid(True)
ax1.set_xlabel('Launch azimuth [deg]')
ax1.set_ylabel('Hyperbolic excess velocity [m/s]')

ax2 = plt.subplot2grid((8,15), (0, 10), rowspan=8, colspan=5, sharey=ax1)
ax2.hist( hev, 40,
          facecolor=colors.cnames['magenta'], edgecolor='black',
          orientation='horizontal', alpha=1.0 )
plt.setp( ax2.get_yticklabels( ), visible=False )
ax2.grid(True)
ax2.set_xlabel('Regolith count')

fig = plt.figure( )
plt.suptitle("Escape characteristic behavior, Launch velocity = 13.0 [m/s]")
data_indices = np.where( initial_velocity == 13.0 )
data_indices = data_indices[ 0 ]
hev = np.sqrt( -muAsteroid / sma[ data_indices ] )

ax1 = plt.subplot2grid((8,15), (0, 0), rowspan=8, colspan=10)
ax1.scatter( azimuth[data_indices], hev )
ax1.grid(True)
ax1.set_xlabel('Launch azimuth [deg]')
ax1.set_ylabel('Hyperbolic excess velocity [m/s]')

ax2 = plt.subplot2grid((8,15), (0, 10), rowspan=8, colspan=5, sharey=ax1)
ax2.hist( hev, 40,
          facecolor=colors.cnames['magenta'], edgecolor='black',
          orientation='horizontal', alpha=1.0 )
plt.setp( ax2.get_yticklabels( ), visible=False )
ax2.grid(True)
ax2.set_xlabel('Regolith count')

## HEV versus inertial launch velocity and launch azimuth angle
fig = plt.figure( )
plt.suptitle("Escape characteristic behavior, Launch velocity = 5.0 [m/s]")
data_indices = np.where( initial_velocity == 5.0 )
data_indices = data_indices[ 0 ]
hev = np.sqrt( -muAsteroid / sma[ data_indices ] )

inertial_vel_x = inertial_vx[data_indices]
inertial_vel_y = inertial_vy[data_indices]
inertial_vel_z = inertial_vz[data_indices]
inertial_vel_mag = np.sqrt( inertial_vel_x**2 + inertial_vel_y**2 + inertial_vel_z**2 )

ax1 = plt.subplot2grid((8,15), (0, 0), rowspan=8, colspan=15)
plot1 = ax1.scatter( azimuth[data_indices], hev, color='black' )
ax1.grid(True)
ax1.set_xlabel('Launch azimuth [deg]')
ax1.set_ylabel('Hyperbolic excess velocity [m/s]')

ax2 = ax1.twinx( )
plotAzimuth = azimuth[data_indices]
plotInertialGuaranteedEscape = inertial_guaranteedEscape[ data_indices ]
plot2, = ax2.plot( plotAzimuth, plotInertialGuaranteedEscape, color='red', linestyle='--', marker='s' )
ax2.set_ylabel('Guaranteed escape speed [m/s]')

ax2.yaxis.label.set_color( plot2.get_color() )

# ax3 = ax1.twiny( )
# ax3.scatter( azimuth[ data_indices ], hev )
# # ax3.set_xticks( ax1.get_xticks( ) )
# # ax3.set_xticklabels( azimuth[data_indices] )
# ax3.set_xlabel('Launch azimuth [deg]')

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
