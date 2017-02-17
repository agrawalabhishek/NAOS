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
Wz = 0.00033118202125129593
mu = 876514

asteroidRotationPeriod = 2.0 * np.pi / Wz
asteroidRotationPeriodHours = asteroidRotationPeriod / ( 60.0 * 60.0 )

## Operations
# Connect to SQLite database.
try:
    database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

phaseAngle = 'N.A.'

## get data for reimpact cases
data1 = pd.read_sql( "SELECT    trajectory_id,                                              \
                                ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth ),                                    \
                                time                                                        \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( crash_flag = 1 );",                                       \
                     database )

data1.columns = [ 'traj_id',                                                                \
                  'vel_mag',                                                                \
                  'launch_azimuth',                                                         \
                  'time' ]

trajectory_id                   = data1[ 'traj_id' ]
velocity_magnitude              = data1[ 'vel_mag' ]
azimuth                         = data1[ 'launch_azimuth' ]
simulation_time                 = data1[ 'time' ]

simulation_time = simulation_time / ( 24.0 * 60.0 * 60.0 )

fig = plt.figure( )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )

## seperate out data for initial launch velocities <= 4
velocityMagnitude_plot = []
azimuth_plot = []
simulationTime_plot = []

for index in range( 0, len( trajectory_id ) ):
    if velocity_magnitude[index] <= 4.0:
        velocityMagnitude_plot.append( velocity_magnitude[index] )
        azimuth_plot.append( azimuth[index] )
        simulationTime_plot.append( simulation_time[index] * 24.0 * 60.0 )

## plot the time to reimpact for velocities <= 4
hexBinPlot = ax1.hexbin( azimuth_plot, simulationTime_plot,                           \
                         C=velocityMagnitude_plot,                                    \
                         cmap='jet', gridsize=100 )
cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax1 )

ax1.grid( True )
ax1.set_xlabel( 'Launch azimuth [deg]' )
ax1.set_ylabel( 'Time to reimpact [hrs]' )
# start, end = ax1.get_ylim( )
# ax1.set_yticks( np.arange( start, end, 10.0 ) )
ax1.set_title( 'Reimpact time versus launch direction and velocity \n Phase angle = '
                + str( phaseAngle ) + ', Asteroid rotation period = '
                + str( asteroidRotationPeriodHours ) + ' [hrs]' )
cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )

## seperate out data for initial launch velocities > 4 and <= 8
velocityMagnitude_plot = []
azimuth_plot = []
simulationTime_plot = []

for index in range( 0, len( trajectory_id ) ):
    if velocity_magnitude[index] <= 8.0 and velocity_magnitude[index] > 4.0:
        velocityMagnitude_plot.append( velocity_magnitude[index] )
        azimuth_plot.append( azimuth[index] )
        simulationTime_plot.append( simulation_time[index] * 24.0 * 60.0 )


## plot the time to reimpact for velocities > 4 and <= 8
hexBinPlot = ax2.hexbin( azimuth_plot, simulationTime_plot,                           \
                         C=velocityMagnitude_plot,                                    \
                         cmap='jet', gridsize=100 )
cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax2 )

ax2.grid( True )
ax2.set_xlabel( 'Launch azimuth [deg]' )
ax2.set_ylabel( 'Time to reimpact [hrs]' )
# start, end = ax2.get_ylim( )
# ax2.set_yticks( np.arange( start, end, 10.0 ) )
ax2.set_title( 'Reimpact time versus launch direction and velocity \n Phase angle = '
                + str( phaseAngle ) )
cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )

## seperate out data for initial launch velocities > 8 and <= 11
velocityMagnitude_plot = []
azimuth_plot = []
simulationTime_plot = []

for index in range( 0, len( trajectory_id ) ):
    if velocity_magnitude[index] <= 11.0 and velocity_magnitude[index] > 8.0:
        velocityMagnitude_plot.append( velocity_magnitude[index] )
        azimuth_plot.append( azimuth[index] )
        simulationTime_plot.append( simulation_time[index] )


## plot the time to reimpact for velocities > 4 and <= 8
hexBinPlot = ax3.hexbin( azimuth_plot, simulationTime_plot,                           \
                         C=velocityMagnitude_plot,                                    \
                         cmap='jet', gridsize=100 )
cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax3 )

ax3.grid( True )
ax3.set_xlabel( 'Launch azimuth [deg]' )
ax3.set_ylabel( 'Time to reimpact [days]' )
# start, end = ax3.get_ylim( )
# ax3.set_yticks( np.arange( start, end, 15.0 ) )
# ax3.set_yscale('log')
ax3.set_title( 'Reimpact time versus launch direction and velocity \n Phase angle = '
                + str( phaseAngle ) )
cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )

## seperate out data for initial launch velocities > 11 and <= 16
velocityMagnitude_plot = []
azimuth_plot = []
simulationTime_plot = []

for index in range( 0, len( trajectory_id ) ):
    if velocity_magnitude[index] <= 16.0 and velocity_magnitude[index] > 11.0:
        velocityMagnitude_plot.append( velocity_magnitude[index] )
        azimuth_plot.append( azimuth[index] )
        simulationTime_plot.append( simulation_time[index] )


## plot the time to reimpact for velocities > 4 and <= 8
hexBinPlot = ax4.hexbin( azimuth_plot, simulationTime_plot,                           \
                         C=velocityMagnitude_plot,                                    \
                         cmap='jet', gridsize=100 )
cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax4 )

ax4.grid( True )
ax4.set_xlabel( 'Launch azimuth [deg]' )
ax4.set_ylabel( 'Time to reimpact [days]' )
start, end = ax4.get_ylim( )
ax4.set_yticks( np.arange( start, end, 10.0 ) )
# ax4.set_yscale('log')
ax4.set_title( 'Reimpact time versus launch direction and velocity \n Phase angle = '
                + str( phaseAngle ) )
cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )

## get data for the escape case
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1, height_ratios = [ 1 ] )
ax3 = plt.subplot( gs[ 0 ] )

data1 = pd.read_sql( "SELECT    trajectory_id,                                              \
                                ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth ),                                    \
                                time                                                        \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 );",                                      \
                     database )

data1.columns = [ 'traj_id',                                                                \
                  'vel_mag',                                                                \
                  'launch_azimuth',                                                         \
                  'time' ]

trajectory_id                   = data1[ 'traj_id' ]
velocity_magnitude              = data1[ 'vel_mag' ]
azimuth                         = data1[ 'launch_azimuth' ]
simulation_time                 = data1[ 'time' ]

simulation_time = simulation_time / ( 24.0 * 60.0 * 60.0 )

hexBinPlot = ax3.hexbin( azimuth, simulation_time,                      \
                         C=velocity_magnitude,                          \
                         cmap='jet', gridsize=100 )
cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax3 )

ax3.grid( True )
ax3.set_xlabel( 'Launch azimuth [deg]' )
ax3.set_ylabel( 'Time to escape [hrs]' )
ax3.set_title( 'Escape time versus launch direction and velocity \n Phase angle = ' + str( phaseAngle ) )
cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )

## show the plot
plt.show( )

# Close SQLite database if it's still open.
if database:
    database.close()

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
