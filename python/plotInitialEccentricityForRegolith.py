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

## function to convert keplerian elements into inertial coordinates
# All angles in the arguments have to be provided in radians (VERY IMPORTANT)
def convertKeplerElementsToCartesianCoordinates( sma,
                                                 eccentricity,
                                                 inclination,
                                                 raan,
                                                 aop,
                                                 ta,
                                                 mu,
                                                 tolerance ):
    # function defination begins here
    if np.abs( eccentricity - 1.0 ) > tolerance:
        semiLatus = sma * ( 1.0 - eccentricity**2 )
    # for parabolic orbits
    else:
        semiLatus = sma

    # compute position and velocity in the perifocal coordinate system
    radius = semiLatus / ( 1.0 + eccentricity * np.cos( ta ) )

    xPositionPerifocal = radius * np.cos( ta )
    yPositionPerifocal = radius * np.sin( ta )

    xVelocityPerifocal = -1.0 * np.sqrt( mu / semiLatus ) * np.sin( ta )
    yVelocityPerifocal = np.sqrt( mu / semiLatus ) * ( eccentricity + np.cos( ta ) )

    # Calculate components of the rotation matrix
    rotation11 = np.cos( raan ) * np.cos( aop ) - np.sin( raan ) * np.sin( aop ) * np.cos( inclination )

    rotation12 = -1.0 * np.cos( raan ) * np.sin( aop ) - np.sin( raan ) * np.cos( aop ) * np.cos( inclination )

    rotation21 = np.sin( raan ) * np.cos( aop ) + np.cos( raan ) * np.sin( aop ) * np.cos( inclination )

    rotation22 = -1.0 * np.sin( raan ) * np.sin( aop ) + np.cos( raan ) * np.cos( aop ) * np.cos( inclination )

    rotation31 = np.sin( aop ) * np.sin( inclination )

    rotation32 = np.cos( aop ) * np.sin( inclination )

    # Compute the cartesian position and velocity in the inertial frame
    xPosition = rotation11 * xPositionPerifocal + rotation12 * yPositionPerifocal

    yPosition = rotation21 * xPositionPerifocal + rotation22 * yPositionPerifocal

    zPosition = rotation31 * xPositionPerifocal + rotation32 * yPositionPerifocal

    xVelocity = rotation11 * xVelocityPerifocal + rotation12 * yVelocityPerifocal

    yVelocity = rotation21 * xVelocityPerifocal + rotation22 * yVelocityPerifocal

    zVelocity = rotation31 * xVelocityPerifocal + rotation32 * yVelocityPerifocal

    # return the values as a tuple
    return ( xPosition, yPosition, zPosition )

## Operations
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593
mu = 876514

## Operations
# Connect to SQLite database.
try:
    # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db" )
    database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

phaseAngle = 'N.A.'

fig = plt.figure( )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )

## get data for the reimpact case
data1 = pd.read_sql( "SELECT    trajectory_id,                                              \
                                ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( crash_flag = 1 );",                                       \
                     database )

data1.columns = [ 'traj_id',                                                                \
                  'vel_mag',                                                                \
                  'launch_azimuth' ]

trajectory_id                   = data1[ 'traj_id' ]
velocity_magnitude              = data1[ 'vel_mag' ]
azimuth                         = data1[ 'launch_azimuth' ]

trajectory_id_list = trajectory_id.tolist( )
trajectory_id_tuple = tuple( trajectory_id_list )

data2 = pd.read_sql( "SELECT    eccentricity,                                                   \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth )                                         \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      trajectory_id IN " + str( trajectory_id_tuple ) + "             \
                     AND        start_flag = 1;",                                               \
                     database )

data2.columns = [ 'eccentricity',                                                       \
                  'vel_mag',                                                            \
                  'launch_azimuth' ]

initial_eccentricity              =  data2[ 'eccentricity' ]
data2_initial_velocity_magnitude  =  data2[ 'vel_mag' ]
data2_initial_azimuth             =  data2[ 'launch_azimuth' ]

## plot data for the reimpact case
unique_reimpact_velocities = np.unique( data2_initial_velocity_magnitude )
gridsize = len( unique_reimpact_velocities )

if gridsize:
    hexBinPlot = ax1.hexbin( data2_initial_azimuth, data2_initial_velocity_magnitude,                 \
                             C=initial_eccentricity,                                                  \
                             reduce_C_function=np.median,                                             \
                             cmap='autumn', gridsize=gridsize, marginals=False )
    cbar = plt.colorbar( hexBinPlot, cmap='autumn', ax=ax1 )

    # ax1.set_yticks( np.arange( 0.0, 14.0, 1.0 ) )
    ax1.set_xticks( np.arange( 0.0, 360.0, 30.0 ) )
    ax1.set_xlim( 0.0, 360.0 )
    ax1.grid( True )
    ax1.set_xlabel( 'Launch azimuth [deg]' )
    ax1.set_ylabel( 'Launch velocity [m/s]' )
    ax1.set_title( 'Reimpact case' )
    cbar.ax.set_ylabel( 'Initial eccentricity' )

    ## scatter plot for the reimpact case
    # get unique velocity values in the escape case
    unique_reimpact_velocities = np.unique( data2_initial_velocity_magnitude )
    colors = plt.cm.Vega20( np.linspace( 0, 1, len( unique_reimpact_velocities ) ) )

    for index in range( 0, len( unique_reimpact_velocities ) ):
        current_velocity = unique_reimpact_velocities[ index ]
        current_velocity_indices = np.where( data2_initial_velocity_magnitude == current_velocity )
        current_velocity_indices = current_velocity_indices[ 0 ]
        plotEccentricity = initial_eccentricity[ current_velocity_indices ]
        plotAzimuth = data2_initial_azimuth[ current_velocity_indices ]
        ax2.scatter( plotAzimuth, plotEccentricity, s=5, c=colors[ index ],                 \
                     edgecolors='face',                                                     \
                     label='$V_{launch}$ = ' + str( current_velocity ) + ' [m/s]' )

    ax2.set_xticks( np.arange( 0.0, 360.0, 30.0 ) )
    ax2.set_xlim( 0.0, 360.0 )
    ax2.grid( True )
    ax2.legend( ).draggable( )
    ax2.set_xlabel( 'Launch azimuth [deg]' )
    ax2.set_ylabel( 'Initial eccentricity' )
    ax2.set_title( 'Reimpact case' )

## get data for the escape case
data1 = pd.read_sql( "SELECT    trajectory_id,                                              \
                                ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 );",                                      \
                     database )

data1.columns = [ 'traj_id',                                                                \
                  'vel_mag',                                                                \
                  'launch_azimuth' ]

trajectory_id                   = data1[ 'traj_id' ]
velocity_magnitude              = data1[ 'vel_mag' ]
azimuth                         = data1[ 'launch_azimuth' ]

trajectory_id_list = trajectory_id.tolist( )
trajectory_id_tuple = tuple( trajectory_id_list )

data2 = pd.read_sql( "SELECT    eccentricity,                                                   \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth )                                         \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      trajectory_id IN " + str( trajectory_id_tuple ) + "             \
                     AND        start_flag = 1;",                                               \
                     database )

data2.columns = [ 'eccentricity',                                                       \
                  'vel_mag',                                                            \
                  'launch_azimuth' ]

initial_eccentricity              =  data2[ 'eccentricity' ]
data2_initial_velocity_magnitude  =  data2[ 'vel_mag' ]
data2_initial_azimuth             =  data2[ 'launch_azimuth' ]

## plot hexbin for escape case
unique_escape_velocities = np.unique( data2_initial_velocity_magnitude )
gridsize = len( unique_escape_velocities )

if gridsize:
    hexBinPlot = ax3.hexbin( data2_initial_azimuth, data2_initial_velocity_magnitude,                 \
                             C=initial_eccentricity,                                                  \
                             reduce_C_function=np.median,                                             \
                             cmap='autumn', gridsize=gridsize, marginals=False )
    cbar = plt.colorbar( hexBinPlot, cmap='autumn', ax=ax3 )

    # ax3.set_yticks( np.arange( 0.0, 14.0, 1.0 ) )
    ax3.set_xticks( np.arange( 0.0, 360.0, 30.0 ) )
    ax3.set_xlim( 0.0, 360.0 )
    ax3.grid( True )
    ax3.set_xlabel( 'Launch azimuth [deg]' )
    ax3.set_ylabel( 'Launch velocity [m/s]' )
    ax3.set_title( 'Escape case' )
    cbar.ax.set_ylabel( 'Initial eccentricity' )

    ## get unique velocity values in the escape case
    unique_escape_velocities = np.unique( data2_initial_velocity_magnitude )
    colors = plt.cm.Vega20( np.linspace( 0, 1, len( unique_escape_velocities ) ) )

    for index in range( 0, len( unique_escape_velocities ) ):
        current_velocity = unique_escape_velocities[ index ]
        current_velocity_indices = np.where( data2_initial_velocity_magnitude == current_velocity )
        current_velocity_indices = current_velocity_indices[ 0 ]
        plotEccentricity = initial_eccentricity[ current_velocity_indices ]
        plotAzimuth = data2_initial_azimuth[ current_velocity_indices ]
        ax4.scatter( plotAzimuth, plotEccentricity, s=5, c=colors[ index ],                 \
                     edgecolors='face',                                                     \
                     label='$V_{launch}$ = ' + str( current_velocity ) + ' [m/s]' )

    ax4.set_xticks( np.arange( 0.0, 360.0, 30.0 ) )
    ax4.set_xlim( 0.0, 360.0 )
    # ax4.set_yticks( np.arange( 12.0, 16.0, 1.0 ) )
    ax4.grid( True )
    ax4.legend( ).draggable( )
    ax4.set_xlabel( 'Launch azimuth [deg]' )
    # ax4.set_ylabel( 'Launch velocity [m/s]' )
    ax4.set_ylabel( 'Initial eccentricity' )
    ax4.set_title( 'Escape case' )
    # cbar.ax.set_ylabel( 'Initial eccentricity' )

##close the database
if database:
    database.close( )

## set global plot title
plt.suptitle( 'Initial eccentricity variation with launch direction and velocity \n Longest edge, \
              Phase angle = ' + str( phaseAngle ) )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## show the plot
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""