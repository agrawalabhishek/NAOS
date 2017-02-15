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
        database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

phaseAngle = 'N.A.'

## get data
data1 = pd.read_sql( "SELECT    trajectory_id,                                              \
                                ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 );",                                      \
                     database )

data1.columns = [ 'traj_id',                                                                \
                  'vel_mag',                                                                \
                  'launch_azimuth' ]

crash_trajectory_id                   = data1[ 'traj_id' ]
crash_velocity_magnitude              = data1[ 'vel_mag' ]
crash_azimuth                         = data1[ 'launch_azimuth' ]

crash_trajectory_id_list = crash_trajectory_id.tolist( )
crash_trajectory_id_tuple = tuple( crash_trajectory_id_list )

data2 = pd.read_sql( "SELECT    sma,                                                            \
                                eccentricity,                                                   \
                                inclination,                                                    \
                                raan,                                                           \
                                aop,                                                            \
                                ta,                                                             \
                                initial_position_x,                                             \
                                initial_position_y,                                             \
                                initial_position_z,                                             \
                                initial_position_magnitude,                                     \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth )                                         \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      trajectory_id IN " + str( crash_trajectory_id_tuple ) + "       \
                     AND        start_flag = 1;",                                               \
                     database )

data2.columns = [ 'sma',                                                                \
                  'eccentricity',                                                       \
                  'inclination',                                                        \
                  'raan',                                                               \
                  'aop',                                                                \
                  'ta',                                                                 \
                  'init_pos_x',                                                         \
                  'init_pos_y',                                                         \
                  'init_pos_z',                                                         \
                  'pos_mag',                                                            \
                  'vel_mag',                                                            \
                  'launch_azimuth' ]

## initial sma and eccentricity for all points that have reimpacted on the surface of the asteroid
crash_initial_sma                       =  data2[ 'sma' ]
crash_initial_eccentricity              =  data2[ 'eccentricity' ]
crash_initial_inclination               =  data2[ 'inclination' ]
crash_initial_raan                      =  data2[ 'raan' ]
crash_initial_aop                       =  data2[ 'aop' ]
crash_initial_ta                        =  data2[ 'ta' ]
launch_initial_position_x               =  data2[ 'init_pos_x' ]
launch_initial_position_y               =  data2[ 'init_pos_y' ]
launch_initial_position_z               =  data2[ 'init_pos_z' ]
data2_crash_initial_position_magnitude  =  data2[ 'pos_mag' ]
data2_crash_initial_velocity_magnitude  =  data2[ 'vel_mag' ]
data2_crash_initial_azimuth             =  data2[ 'launch_azimuth' ]

## convert the angles from degrees to radians
crash_initial_inclination   = crash_initial_inclination * np.pi / 180.0
crash_initial_raan          = crash_initial_raan * np.pi / 180.0
crash_initial_aop           = crash_initial_aop * np.pi / 180.0

initialPericentres = crash_initial_sma * ( 1.0 - crash_initial_eccentricity )

## make the true anomaly zero to get the pericentre points
crash_initial_ta = np.zeros( len( crash_initial_sma ) )

## convert orbital elements at true-anomaly=0 to cartesian coordinates
# this would give the location of the initial periapsis in the x-y-z coordinate frame
# note that at time t=0 the inertial and body frame are alligned so the position vector would be the
# same for both frames

# Compute semi-latus rectum
# semiLatus for non parabolic orbits (eccentricity != 1)
tolerance = 1.0e-15
# semiLatus = []
# radius = []
# xPositionPerifocal = []
# yPositionPerifocal = []
# xVelocityPerifocal = []
# yVelocityPerifocal = []
# rotation11 = []
# rotation12 = []
# rotation21 = []
# rotation22 = []
# rotation31 = []
# rotation32 = []
xPosition                   = []
yPosition                   = []
zPosition                   = []
checker                     = []
checkerPlot                 = []
crashAzimuthPlot            = []
crashPericentrePlot         = []
crashInitialVelocityPlot    = []

sma_test             = 6787746.891;
eccentricity_test    = 0.000731104;
inclination_test     = 51.68714486 * np.pi / 180.0;
raan_test            = 127.5486706 * np.pi / 180.0;
aop_test             = 74.21987137 * np.pi / 180.0;
ta_test              = 24.10027677 * np.pi / 180.0;
gravParameter_test   = 3.98600441e14;

( xPosition_test, yPosition_test, zPosition_test )= convertKeplerElementsToCartesianCoordinates(
                                                                               sma_test,
                                                                               eccentricity_test,
                                                                               inclination_test,
                                                                               raan_test,
                                                                               aop_test,
                                                                               ta_test,
                                                                               gravParameter_test,
                                                                               tolerance )

print "x position test value = " + str( xPosition_test )
print "y position test value = " + str( yPosition_test )
print "z position test value = " + str( zPosition_test )

for index in range( 0, len( crash_initial_sma ) ):
    ( xCoordinate, yCoordinate, zCoordinate ) = convertKeplerElementsToCartesianCoordinates(
                                                       crash_initial_sma[index],
                                                       crash_initial_eccentricity[index],
                                                       crash_initial_inclination[index],
                                                       crash_initial_raan[index],
                                                       crash_initial_aop[index],
                                                       crash_initial_ta[index],
                                                       mu,
                                                       tolerance )

    xPosition.append( xCoordinate )
    yPosition.append( yCoordinate )
    zPosition.append( zCoordinate )

    checker.append( ( xPosition[index]**2 / alpha**2 ) + ( yPosition[index]**2 / beta**2 ) + ( zPosition[index]**2 / gamma**2 ) )
    if checker[index] > 1.0:
        checkerPlot.append( checker[index] )
        crashAzimuthPlot.append( data2_crash_initial_azimuth[index] )
        crashPericentrePlot.append( initialPericentres[index] )
        crashInitialVelocityPlot.append( data2_crash_initial_velocity_magnitude[index] )

print checkerPlot

if checkerPlot:
    ## plot the initial pericentre values for the reimpacted particles
    fig = plt.figure( )
    gs = gridspec.GridSpec( 3, 1, height_ratios = [ 1, 1, 1 ] )
    ax1 = plt.subplot( gs[ 0 ] )
    # ax2 = plt.subplot( gs[ 1 ], projection='3d' )
    ax2 = plt.subplot( gs[ 1 ] )
    ax3 = plt.subplot( gs[ 2 ] )

    hexBinPlot = ax1.hexbin( data2_crash_initial_azimuth, initialPericentres, data2_crash_initial_velocity_magnitude, cmap='jet' )
    cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax1 )

    ax1.hlines( data2_crash_initial_position_magnitude[ 0 ], 0.0, 359.0, colors=colors.cnames['purple'] )

    ax1.set_ylim( 0, 8500 )
    ax1.grid( True )
    ax1.set_xlabel( 'Launch azimuth [deg]' )
    ax1.set_ylabel( 'Initial pericentre distance [m]' )
    ax1.set_title( 'Initial Pericentre variation with launch direction and velocity \n Phase angle = ' + str( phaseAngle ) + ', Escape case' )
    cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )

    alpha = 20000.0
    beta = 7000.0
    gamma = 7000.0

    hexBinPlot2 = ax2.hexbin( crashAzimuthPlot, crashPericentrePlot, crashInitialVelocityPlot, cmap='inferno' )
    cbar = plt.colorbar( hexBinPlot2, cmap='inferno', ax=ax2 )
    ax2.set_xlim( 0, 360 )
    ax2.set_ylim( 0, 8500 )
    ax2.set_xlabel( 'Launch azimuth [deg]' )
    ax2.set_ylabel( 'Initial pericentre distance [m]' )
    cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )
    ax2.set_title( 'Initial pericentre points outside the ellipsoid' )
    ax2.grid( True )

    ## plot the values of the solutions to the ellipsoid equation
    # hexBinPlot3 = ax3.hexbin( crashAzimuthPlot, checkerPlot, crashInitialVelocityPlot, cmap='inferno' )
    hexBinPlot3 = ax3.hexbin( crashAzimuthPlot, crashPericentrePlot, checkerPlot, cmap='inferno' )
    cbar = plt.colorbar( hexBinPlot3, cmap='inferno', ax=ax3 )
    # ax3.set_xlim( 0, 360 )
    # ax3.set_ylim( 0, 8500 )
    ax3.set_xlabel( 'Launch azimuth [deg]' )
    ax3.set_ylabel( 'Initial pericentre distance [m]' )
    cbar.ax.set_ylabel( 'Ellipsoid equation solution' )
    ax3.set_title( 'Initial pericentre points outside the ellipsoid' )
    ax3.grid( True )
else:
    ## plot the initial pericentre values for the reimpacted particles
    fig = plt.figure( )
    gs = gridspec.GridSpec( 2, 1, height_ratios = [ 1, 1 ] )
    ax1 = plt.subplot( gs[ 0 ] )
    ax2 = plt.subplot( gs[ 1 ], projection = '3d' )

    hexBinPlot = ax1.hexbin( data2_crash_initial_azimuth, initialPericentres,                   \
                             C=data2_crash_initial_velocity_magnitude,                          \
                             cmap='jet', gridsize=50 )
    cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax1 )

    ax1.hlines( data2_crash_initial_position_magnitude[ 0 ], 0.0, 359.0, colors=colors.cnames['purple'] )

    ax1.set_ylim( 0, 8500 )
    ax1.grid( True )
    ax1.set_xlabel( 'Launch azimuth [deg]' )
    ax1.set_ylabel( 'Initial pericentre distance [m]' )
    ax1.set_title( 'Initial Pericentre variation with launch direction and velocity \n Phase angle = ' + str( phaseAngle ) + ', Escape case' )
    cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )

    # plot the 3d ellipsoid wireframe and siplay pericentres in it
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
    ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
    ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

    surf = ax2.plot_wireframe( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                               rstride=5, cstride=5 )

    ax2.scatter( launch_initial_position_x, launch_initial_position_y, launch_initial_position_z,
                 zdir = 'z', c = colors.cnames['green'], label = 'launch site' )

    ax2.scatter( xPosition, yPosition, zPosition,
                 zdir = 'z', c = colors.cnames['red'], label = 'initial pericentres' )

    ax2.set_xlabel( 'x [m]' )
    ax2.set_ylabel( 'y [m]' )
    ax2.set_zlabel( 'z [m]' )
    ax2.set_title( 'Initial pericentre locations' )

## show the plot
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
