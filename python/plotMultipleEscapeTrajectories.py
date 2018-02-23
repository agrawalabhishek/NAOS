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

muAsteroid = 876514

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
commonPath = "/media/abhishek/Ashish/Thesis_Simulation_Databases/trailing_edge_perturbations_CDE/"
try:
        # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
        #                             + "multiple_launch_velocity/"
        #                             + "simulation_time_9_months/"
        #                             + "longestEdge.db")
        # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
        #                            + "multiple_launch_velocity_with_perturbations/"
        #                            + "simulation_time_9_months/"
        #                            + "3.2Density_1cmRadius/longestEdgePerturbations.db")
        database = sqlite3.connect("../data/regolith_launched_from_trailing_edge/"
                                   + "multiple_launch_velocity_with_perturbations/"
                                   + "simulation_time_9_months/"
                                   + "trailingEdge_3P2Density_1cmRadius.db")
        # database = sqlite3.connect(commonPath + "7.5Density_5cmRadius/trailingEdge_7P5Density_5cmRadius.db")
        # database = sqlite3.connect(commonPath + "3.2Density_1cmRadius/trailingEdge_3P2Density_1cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

azimuthRange = np.arange( 270.0, 285.0, 1.0 )
# appendAzimuth = np.arange( 0.0, 10.0, 1.0 )
# azimuthRange = np.append( azimuthRange, appendAzimuth )
azimuthRange = tuple( azimuthRange.tolist( ) )
# azimuthRange = ( 270.0 )

print "Extracting data now...\n"

# data = pd.read_sql( "SELECT     trajectory_id,                                              \
#                                 ROUND( launch_azimuth ),                                    \
#                                 sma,                                                        \
#                                 time                                                        \
#                      FROM       regolith_trajectory_results                                 \
#                      WHERE      ROUND( escape_flag ) = 1                                    \
#                      AND        ROUND( initial_velocity_magnitude ) = 9.0                   \
#                      AND        ROUND( launch_azimuth ) IN " + str( azimuthRange ) + ";",   \
#                      database )

data = pd.read_sql( "SELECT     trajectory_id,                                              \
                                ROUND( launch_azimuth ),                                    \
                                sma,                                                        \
                                time                                                        \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ROUND( escape_flag ) = 1                                    \
                     AND        ROUND( initial_velocity_magnitude ) = 9.0                   \
                     AND        ROUND( initial_solar_phase_angle ) = 225.0                  \
                     AND        ROUND( launch_azimuth ) IN " + str( azimuthRange ) + ";",   \
                     database )

data.columns = [ 'trajectory_id', 'launch_azimuth', 'sma', 'time' ]

trajectory_id = data[ 'trajectory_id' ]
sma           = data[ 'sma' ]
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
                                sma,                                                        \
                                eccentricity,                                               \
                                total_energy,                                               \
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
                  'sma',                                                 \
                  'eccentricity',                                        \
                  'total_energy',                                        \
                  'launch_azimuth' ]

x                   = data1[ 'x' ]
y                   = data1[ 'y' ]
z                   = data1[ 'z' ]
bodyFrame_x         = data1[ 'bodyFrame_x' ]
bodyFrame_y         = data1[ 'bodyFrame_y' ]
bodyFrame_z         = data1[ 'bodyFrame_z' ]
initial_velocity    = data1[ 'init_vel_mag' ]
t                   = data1[ 'time' ]
fullSMA             = data1[ 'sma' ]
fullEcc             = data1[ 'eccentricity' ]
fullEnergy          = data1[ 'total_energy' ]
azimuth             = data1[ 'launch_azimuth' ]

xPositionStart      = data1[ 'init_pos_x' ]
yPositionStart      = data1[ 'init_pos_y' ]
zPositionStart      = data1[ 'init_pos_z' ]

initial_velocity = initial_velocity.tolist( )

if database:
    database.close( )

# Extraction time
end_time = time.time( )
print "Data extraction time: " + str("{:,g}".format(end_time - start_time)) + "s\n"

print "Plotting data now...\n"

viewAngles = [ 45 ] # for the mayavi 3D view

## HEV and escape count
fig = plt.figure( )
hev = np.sqrt( -muAsteroid / sma )

ax1 = plt.subplot2grid((8,10), (0, 0), rowspan=8, colspan=10)
ax1.scatter( launchAzim, hev )
ax1.grid( True )
ax1.set_xlabel( 'Launch azimuth [deg]' )
ax1.set_ylabel( 'Hyperbolic excess velocity [m/s]' )

## Plot all escape trajectories for the given azimuth angle range and launch velocity
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

plotEllipse( alpha, beta, ax1 )
uniqueAzimuths = np.unique( azimuth )
for azimuthIndex in range( 0, len( uniqueAzimuths ) ):
    currentAzimuth = uniqueAzimuths[ azimuthIndex ]
    data_indices = np.where( azimuth == currentAzimuth )
    data_indices = data_indices[ 0 ]
    trajectory_x = x[data_indices]
    trajectory_y = y[data_indices]
    trajectory_z = z[data_indices]

    ax1.plot( trajectory_x, trajectory_y,
              linewidth=1.0,
              label="Azimuth = " + str( currentAzimuth ) + " [deg]" )

ax1.grid(True)
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax1.legend( ).draggable( )

## SMA and eccentricity
fig = plt.figure( figsize=(8, 8) )
gs = gridspec.GridSpec( 2, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
plt.suptitle("Escape characteristic behavior, Launch velocity = " + str(initial_velocity[ 0 ]) + " [m/s]")

uniqueAzimuths = np.unique( azimuth )
for azimuthIndex in range( 0, len( uniqueAzimuths ) ):
    currentAzimuth = uniqueAzimuths[ azimuthIndex ]
    data_indices = np.where( azimuth == currentAzimuth )
    data_indices = data_indices[ 0 ]
    trajectory_x = x[data_indices]
    trajectory_y = y[data_indices]
    trajectory_z = z[data_indices]

    plot_sma = fullSMA[ data_indices ]
    plot_ecc = fullEcc[ data_indices ]
    plot_energy = fullEnergy[ data_indices ]
    plotEnergy = plot_energy.tolist( )
    plot_time = t[ data_indices ]

    ax1.plot( plot_time, plot_ecc, lw=1.0,
              label=str( currentAzimuth ) )

    ax2.plot( plot_time, plot_energy, lw=1.0,
              label=str( currentAzimuth ) )

ax1.grid(True)
ax1.set_xlabel('Time [s]')
ax1.set_xscale('log')
# ax1.set_ylabel('Initial Energy [$m^2/s^2$]')
ax1.set_ylabel('Eccentricity')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax1.legend( title='Launch Azimuth [deg]' ).draggable( )

ax2.grid(True)
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Energy [$m^2/s^2$]')
# ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax2.legend( title='Launch Azimuth [deg]', markerscale=10.0 ).draggable( )

plt.show( )
sys.exit( )

## SMA and eccentricity
fig = plt.figure( figsize=(8, 8) )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )
plt.suptitle("Escape characteristic behavior, Launch velocity = " + str(initial_velocity[ 0 ]) + " [m/s]")

uniqueAzimuths = np.unique( azimuth )
for azimuthIndex in range( 0, len( uniqueAzimuths ) ):
    currentAzimuth = uniqueAzimuths[ azimuthIndex ]
    data_indices = np.where( azimuth == currentAzimuth )
    data_indices = data_indices[ 0 ]
    trajectory_x = x[data_indices]
    trajectory_y = y[data_indices]
    trajectory_z = z[data_indices]

    plot_sma = fullSMA[ data_indices ]
    plot_ecc = fullEcc[ data_indices ]
    plot_energy = fullEnergy[ data_indices ]
    plotEnergy = plot_energy.tolist( )
    plot_time = t[ data_indices ]

    if currentAzimuth <= 9.0:
        # ax1.scatter( currentAzimuth, plotEnergy[ 0 ] )
        ax1.plot( plot_time, plot_ecc, lw=1.0,
                  label=str( currentAzimuth ) )

        ax3.plot( plot_time, plot_energy, lw=1.0,
                  label=str( currentAzimuth ) )
    else:
        # ax1.scatter( currentAzimuth, plotEnergy[ 0 ] )
        ax2.plot( plot_time, plot_ecc, lw=1.0,
                  label=str( currentAzimuth ) )

        ax4.plot( plot_time, plot_energy, lw=1.0,
                  label=str( currentAzimuth ) )

ax1.grid(True)
ax1.set_xlabel('Time [s]')
ax1.set_xscale('log')
# ax1.set_ylabel('Initial Energy [$m^2/s^2$]')
ax1.set_ylabel('Eccentricity')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax1.legend( title='Launch Azimuth [deg]' ).draggable( )

ax3.grid(True)
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Energy [$m^2/s^2$]')
# ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax3.legend( title='Launch Azimuth [deg]', markerscale=10.0 ).draggable( )

ax2.grid(True)
ax2.set_xlabel('Time [s]')
ax2.set_xscale('log')
# ax2.set_ylabel('Initial Energy [$m^2/s^2$]')
ax2.set_ylabel('Eccentricity')
# ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax2.legend( title='Launch Azimuth [deg]' ).draggable( )

ax4.grid(True)
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_xlabel('Time [s]')
ax4.set_ylabel('Energy [$m^2/s^2$]')
# ax4.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useoffset=False)
ax4.legend( title='Launch Azimuth [deg]', markerscale=10.0 ).draggable( )

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
