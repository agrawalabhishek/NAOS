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
                     WHERE      ( escape_flag = 1 );",                                       \
                     database )

data1.columns = [ 'traj_id',                                                                \
                  'vel_mag',                                                                \
                  'launch_azimuth' ]

trajectory_id                   = data1[ 'traj_id' ]
velocity_magnitude              = data1[ 'vel_mag' ]
azimuth                         = data1[ 'launch_azimuth' ]

trajectory_id_list = trajectory_id.tolist( )
trajectory_id_tuple = tuple( trajectory_id_list )

data2 = pd.read_sql( "SELECT    sma,                                                            \
                                eccentricity,                                                   \
                                inclination,                                                    \
                                raan,                                                           \
                                aop,                                                            \
                                ta,                                                             \
                                total_energy,                                                   \
                                initial_position_x,                                             \
                                initial_position_y,                                             \
                                initial_position_z,                                             \
                                initial_position_magnitude,                                     \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth )                                         \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      trajectory_id IN " + str( trajectory_id_tuple ) + "             \
                     AND        start_flag = 1;",                                               \
                     database )

data2.columns = [ 'sma',                                                                \
                  'eccentricity',                                                       \
                  'inclination',                                                        \
                  'raan',                                                               \
                  'aop',                                                                \
                  'ta',                                                                 \
                  'total_energy',                                                       \
                  'init_pos_x',                                                         \
                  'init_pos_y',                                                         \
                  'init_pos_z',                                                         \
                  'pos_mag',                                                            \
                  'vel_mag',                                                            \
                  'launch_azimuth' ]

## initial sma and eccentricity for all points that have reimpacted on the surface of the asteroid
initial_sma                       =  data2[ 'sma' ]
initial_eccentricity              =  data2[ 'eccentricity' ]
initial_inclination               =  data2[ 'inclination' ]
initial_raan                      =  data2[ 'raan' ]
initial_aop                       =  data2[ 'aop' ]
initial_ta                        =  data2[ 'ta' ]
initial_energy                    =  data2[ 'total_energy' ]
launch_initial_position_x         =  data2[ 'init_pos_x' ]
launch_initial_position_y         =  data2[ 'init_pos_y' ]
launch_initial_position_z         =  data2[ 'init_pos_z' ]
data2_initial_position_magnitude  =  data2[ 'pos_mag' ]
data2_initial_velocity_magnitude  =  data2[ 'vel_mag' ]
data2_initial_azimuth             =  data2[ 'launch_azimuth' ]

## convert the angles from degrees to radians
initial_inclination   = initial_inclination * np.pi / 180.0
initial_raan          = initial_raan * np.pi / 180.0
initial_aop           = initial_aop * np.pi / 180.0

## plot with initial velocity on the colorbar
fig = plt.figure( )
gs = gridspec.GridSpec( 2, 1, height_ratios = [ 1, 1 ] )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )

hexBinPlot = ax1.hexbin( initial_eccentricity, initial_energy,                        \
                         C=data2_initial_velocity_magnitude,                          \
                         cmap='jet', gridsize=100 )
cbar = plt.colorbar( hexBinPlot, cmap='jet', ax=ax1 )

ax1.grid( True )
ax1.set_xlabel( 'Eccentricity' )
ax1.set_ylabel( 'Initial energy' + '$[m^2/s^2]$' )
ax1.set_title( 'Initial eccentricity and energy variation with launch velocity \n Phase angle = ' + str( phaseAngle ) + ', Escape case' )
cbar.ax.set_ylabel( 'Regolith launch velocity [m/s]' )

## plot with initial launch direction on the colorbar
hexBinPlot2 = ax2.hexbin( initial_eccentricity, initial_energy,                       \
                         C=data2_initial_azimuth,                                     \
                         cmap='jet', gridsize=100 )
cbar = plt.colorbar( hexBinPlot2, cmap='jet', ax=ax2 )

ax2.grid( True )
ax2.set_xlabel( 'Eccentricity' )
ax2.set_ylabel( 'Initial energy' + '$[m^2/s^2]$' )
ax2.set_title( 'Initial eccentricity and energy variation with launch direction \n Phase angle = ' + str( phaseAngle ) + ', Escape case' )
cbar.ax.set_ylabel( 'Regolith launch azimuth [deg]' )

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
