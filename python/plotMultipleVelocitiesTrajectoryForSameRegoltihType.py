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

def extractDataFromSQL( databaseConnect ):
    # general function for data extraction
    data = pd.read_sql( "SELECT    inertial_position_x,                                         \
                                   inertial_position_y,                                         \
                                   inertial_position_z,                                         \
                                   position_x,                                                  \
                                   position_y,                                                  \
                                   position_z,                                                  \
                                   ROUND( initial_velocity_magnitude ),                         \
                                   ROUND( launch_azimuth ),                                     \
                                   ROUND( initial_solar_phase_angle )                           \
                        FROM       regolith_trajectory_results                                  \
                        WHERE      ROUND( initial_velocity_magnitude ) IN (8.0, 9.0, 10.0)      \
                        AND        ROUND( initial_solar_phase_angle ) = 45.0                    \
                        AND        ROUND( launch_azimuth ) = 270.0;",                           \
                        databaseConnect )

    data.columns = [ 'inertial_pos_x',                                                          \
                     'inertial_pos_y',                                                          \
                     'inertial_pos_z',                                                          \
                     'x',                                                                       \
                     'y',                                                                       \
                     'z',                                                                       \
                     'init_vel_mag',                                                            \
                     'launch_azimuth',                                                          \
                     'init_solar_phase' ]

    x_data                   = data[ 'x' ]
    y_data                   = data[ 'y' ]
    z_data                   = data[ 'z' ]
    initial_velocity_data    = data[ 'init_vel_mag' ]
    azimuth_data             = data[ 'launch_azimuth' ]
    init_solar_phase         = data[ 'init_solar_phase' ]

    xPositionInertial        = data[ 'inertial_pos_x' ]
    yPositionInertial        = data[ 'inertial_pos_y' ]
    zPositionInertial        = data[ 'inertial_pos_z' ]

    return x_data, y_data, z_data, initial_velocity_data, azimuth_data, init_solar_phase,       \
           xPositionInertial, yPositionInertial, zPositionInertial

def areaToMassRatio( radius, density ):
    crossSectionalArea = np.pi * radius**2
    mass = ( 4.0 / 3.0 ) * np.pi * radius**3 * density
    ratio = crossSectionalArea / mass
    return ratio

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

## Constants for the asteroid
alpha = 20000.0
beta = 7000.0
gamma = 7000.0

commonPath = "/media/abhishek/Ashish/Thesis_Simulation_Databases/trailing_edge_perturbations_CDE/"
commonPath2 = "../data/regolith_launched_from_trailing_edge/multiple_launch_velocity_with_perturbations/simulation_time_9_months/"

try:
    database1 = sqlite3.connect( commonPath + "3.2Density_1cmRadius/trailingEdge_3P2Density_1cmRadius.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

# try:
#     database2 = sqlite3.connect( commonPath + "7.5Density_1cmRadius/trailingEdge_7P5Density_1cmRadius.db")

# except sqlite3.Error, e:
#         print "Error %s:" % e.args[0]
#         sys.exit(1)

# try:
#     database3 = sqlite3.connect( commonPath + "3.2Density_5cmRadius/trailingEdge_3P2Density_5cmRadius.db")

# except sqlite3.Error, e:
#         print "Error %s:" % e.args[0]
#         sys.exit(1)

# try:
#     database4 = sqlite3.connect( commonPath + "7.5Density_5cmRadius/trailingEdge_7P5Density_5cmRadius.db")

# except sqlite3.Error, e:
#         print "Error %s:" % e.args[0]
#         sys.exit(1)

fig = plt.figure( figsize=( 7, 7 ) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

ratio = areaToMassRatio( 1.0e-2, 3.2e3 )
label1 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'
( x_data1, y_data1, z_data1, initial_velocity_data1, azimuth_data1, init_solar_phase_data1,
    xPositionInertial_data1, yPositionInertial_data1, zPositionInertial_data1 )                 \
        = extractDataFromSQL( database1 )

# ratio = areaToMassRatio( 1.0e-2, 7.5e3 )
# label2 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'
# ( x_data2, y_data2, z_data2, initial_velocity_data2, azimuth_data2, init_solar_phase_data2,
#     xPositionInertial_data2, yPositionInertial_data2, zPositionInertial_data2 )                 \
#         = extractDataFromSQL( database2 )

# ratio = areaToMassRatio( 5.0e-2, 3.2e3 )
# label3 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'
# ( x_data3, y_data3, z_data3, initial_velocity_data3, azimuth_data3, init_solar_phase_data3,
#     xPositionInertial_data3, yPositionInertial_data3, zPositionInertial_data3 )                 \
#         = extractDataFromSQL( database3 )

# ratio = areaToMassRatio( 5.0e-2, 7.5e3 )
# label4 = 'A/M = ' + str(ratio) + ' [$m^2/kg$]'
# ( x_data4, y_data4, z_data4, initial_velocity_data4, azimuth_data4, init_solar_phase_data4,
#     xPositionInertial_data4, yPositionInertial_data4, zPositionInertial_data4 )                 \
#         = extractDataFromSQL( database4 )

plotEllipse( alpha, beta, ax1 )

uniqueVelocities = np.unique( initial_velocity_data1 )

for index in range(0, len(uniqueVelocities)):
    currentVel = uniqueVelocities[ index ]
    dataIndices = np.where( initial_velocity_data1 == currentVel )
    dataIndices = dataIndices[ 0 ]

    ax1.plot( xPositionInertial_data1[dataIndices], yPositionInertial_data1[dataIndices],
                lw=1.0, label=str(currentVel) )
    # ax1.plot( xPositionInertial_data2, yPositionInertial_data2, lw=1.0, label=label2 )
    # ax1.plot( xPositionInertial_data3, yPositionInertial_data3, lw=1.0, label=label3 )
    # ax1.plot( xPositionInertial_data4, yPositionInertial_data4, lw=1.0, label=label4 )

ax1.text( xPositionInertial_data1[0], yPositionInertial_data1[0],
            'start', fontsize=10, color='black' )

ax1.grid(True)
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_title(label1 + ", Azimuth = " + str(azimuth_data1[ 0 ]) + " [deg], "
                + "Solar phase = " + str(init_solar_phase_data1[ 0 ]) + " [deg]" )

ax1.legend( title="Launch velocity [m/s]").draggable( )
ax1.axis('equal')

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
