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
from decimal import Decimal

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
from matplotlib import animation
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

# asteroid's characteristic constants
alpha = 20000.0
beta = 7000.0
gamma = 7000.0
Wz = 0.00033118202125129593
mu = 876514

# Set up the figure
fig = plt.figure( figsize=(10, 10) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

# Connect to SQLite database.
savePath1 = "../data/regolith_launched_from_longest_edge/"
savePath2 = "multiple_launch_velocity_with_perturbations/"
savePath3 = "simulation_time_9_months/"
savePath4 = "3.2Density_1cmSize/"
savePath5 = "test.mp4"
savePath = savePath1 + savePath2 + savePath3 + savePath4 + savePath5

try:
    database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                               + "multiple_launch_velocity_with_perturbations/"
                               + "simulation_time_9_months/"
                               + "3.2Density_1cmSize/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

lowerTimeDays = 0.0
lowerTime = lowerTimeDays * 24.0 * 60.0 * 60.0

upperTimeDays = 270.0
upperTime = upperTimeDays * 24.0 * 60.0 * 60.0

data = pd.read_sql( "SELECT     ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth ),                                        \
                                ROUND( initial_solar_phase_angle ),                             \
                                solar_phase_angle,                                              \
                                position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                inertial_position_x,                                            \
                                inertial_position_y,                                            \
                                inertial_position_z,                                            \
                                inertial_velocity_x,                                            \
                                inertial_velocity_y,                                            \
                                inertial_velocity_z,                                            \
                                srp_x,                                                          \
                                srp_y,                                                          \
                                srp_z,                                                          \
                                solarTide_x,                                                    \
                                solarTide_y,                                                    \
                                solarTide_z,                                                    \
                                gravAcc_x,                                                      \
                                gravAcc_y,                                                      \
                                gravAcc_z,                                                      \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( launch_azimuth ) = 45.0                                  \
                     AND        ROUND( initial_velocity_magnitude ) = 10.0                      \
                     AND        ROUND( initial_solar_phase_angle ) = 315.0                      \
                     AND        time >= " + str(lowerTime)
                                + " AND time <= " + str(upperTime) + " ;",                      \
                     database )

if database:
    database.close( )

data.columns = [ 'initial_velocity_magnitude',                                                  \
                 'launch_azimuth',                                                              \
                 'initial_solar_phase_angle',                                                   \
                 'solar_phase_angle',                                                           \
                 'position_x',                                                                  \
                 'position_y',                                                                  \
                 'position_z',                                                                  \
                 'inertial_position_x',                                                         \
                 'inertial_position_y',                                                         \
                 'inertial_position_z',                                                         \
                 'inertial_velocity_x',                                                         \
                 'inertial_velocity_y',                                                         \
                 'inertial_velocity_z',                                                         \
                 'srp_x',                                                                       \
                 'srp_y',                                                                       \
                 'srp_z',                                                                       \
                 'solarTide_x',                                                                 \
                 'solarTide_y',                                                                 \
                 'solarTide_z',                                                                 \
                 'gravAcc_x',                                                                   \
                 'gravAcc_y',                                                                   \
                 'gravAcc_z',                                                                   \
                 'time' ]

initialVelocityMagnitude    = data[ 'initial_velocity_magnitude' ]
launchAzimuth               = data[ 'launch_azimuth' ]
solarPhase                  = data[ 'initial_solar_phase_angle' ]
changingSolarPhase          = data[ 'solar_phase_angle' ]
position_x                  = data[ 'position_x' ]
position_y                  = data[ 'position_y' ]
position_z                  = data[ 'position_z' ]
inertial_position_x         = data[ 'inertial_position_x' ]
inertial_position_y         = data[ 'inertial_position_y' ]
inertial_position_z         = data[ 'inertial_position_z' ]
inertial_velocity_x         = data[ 'inertial_velocity_x' ]
inertial_velocity_y         = data[ 'inertial_velocity_y' ]
inertial_velocity_z         = data[ 'inertial_velocity_z' ]
srp_x                       = data[ 'srp_x' ]
srp_y                       = data[ 'srp_y' ]
srp_z                       = data[ 'srp_z' ]
solarTide_x                 = data[ 'solarTide_x' ]
solarTide_y                 = data[ 'solarTide_y' ]
solarTide_z                 = data[ 'solarTide_z' ]
gravAcc_x                   = data[ 'gravAcc_x' ]
gravAcc_y                   = data[ 'gravAcc_y' ]
gravAcc_z                   = data[ 'gravAcc_z' ]
t                           = data[ 'time' ]

print "Processing data now...\n"

print "Solar phase at start of the segment = " + str( changingSolarPhase[ 0 ] ) + " [deg]"

plt.suptitle( '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]', fontsize=12 )

# Draw the entire trajectory path
ax1.plot( inertial_position_x, inertial_position_z, linewidth=1, color=colors.cnames['purple'] )

def init( ):
    ax1.plot( inertial_position_x[ 0 ], inertial_position_z[ 0 ] )

def animate( i ):
    timeText.set_text( time_template % ( timeData[ i ] ) )
    velocityText.set_text( velocity_template % ( velocityData[ i ] ) )
    rangeText.set_text( range_template % ( rangeData[ i ] ) )
    srpText.set_text( srp_template % ( srpData_SciNot[ i ] ) )
    solarTideText.set_text( solarTide_template % ( solarTideData_SciNot[ i ] ) )
    gravAccText.set_text( gravAcc_template % ( gravAcc_SciNot[ i ] ) )
    point.set_data( coordinate[ 0 ][ i ], coordinate[ 1 ][ i ] )
    return point, timeText, velocityText, rangeText, srpText, solarTideText, gravAccText

### main segment of the code
## data indices
steps = 100
data_indices = range( 0, len( inertial_position_x ), steps )

## define the text showing change in time value
time_template = 'Time = %.2f [days]'
# define text placement coordinates
# xTop_text = max( inertial_position_x ) + 0.25e5
# yTop_text = max( inertial_position_z )
xTop_text = -0.5e6
zTop_text = 11.0e5
timeText = ax1.text( xTop_text, zTop_text, '', fontsize=12 )
timeData = t[ data_indices ] / ( 24.0 * 60.0 * 60.0 )
timeData = timeData.tolist( )

## define the text handle showing velocity magnitude of the regolith
velocity_template = 'Velocity = %0.3f [m/s]'
# define velocity placement coordinates
xTop_velocity = xTop_text
zTop_velocity = zTop_text - 0.50e5
velocityText = ax1.text( xTop_velocity, zTop_velocity, '', fontsize=12 )
Vx = inertial_velocity_x[data_indices]
Vy = inertial_velocity_y[data_indices]
Vz = inertial_velocity_z[data_indices]
velocityData = np.sqrt( Vx**2 + Vy**2 + Vz**2 )
velocityData = velocityData.tolist( )

## define the text handle showing range magnitude of the regolith
range_template = 'Range = %0.3f [m]'
# define velocity placement coordinates
xTop_range = xTop_velocity
zTop_range = zTop_velocity - 0.50e5
rangeText = ax1.text( xTop_range, zTop_range, '', fontsize=12 )
rangeData = np.sqrt( inertial_position_x[data_indices]**2 + inertial_position_y[data_indices]**2 + inertial_position_z[data_indices]**2 )
rangeData = rangeData.tolist( )

## define the text handle showing various accelerations acting on the regolith
# SRP
srp_template = 'SRP = %s $[m/s^2]$'
xTop_srp = xTop_range
zTop_srp = zTop_range - 0.50e5
srpText = ax1.text( xTop_srp, zTop_srp, '', fontsize=12 )
xSRP = srp_x[data_indices]
ySRP = srp_y[data_indices]
zSRP = srp_z[data_indices]
srpData = np.sqrt( xSRP**2 + ySRP**2 + zSRP**2 )
srpData = srpData.tolist( )
srpData_SciNot = []
for index in range( 0, len( srpData ) ):
    srpData_SciNot.append( "{:.2e}".format( srpData[index] ) )

# Solar tide
solarTide_template = 'Solar tide = %s $[m/s^2]$'
xTop_solarTide = xTop_srp
zTop_solarTide = zTop_srp - 0.50e5
solarTideText = ax1.text( xTop_solarTide, zTop_solarTide, '', fontsize=12 )
xSTBE = solarTide_x[data_indices]
ySTBE = solarTide_y[data_indices]
zSTBE = solarTide_z[data_indices]
solarTideData = np.sqrt( xSTBE**2 + ySTBE**2 + zSTBE**2 )
solarTideData = solarTideData.tolist( )
solarTideData_SciNot = []
for index in range( 0, len( solarTideData ) ):
    solarTideData_SciNot.append( "{:.2e}".format( solarTideData[index] ) )

# Gravity
gravAcc_template = 'Grav. = %s $[m/s^2]$'
xTop_gravAcc = xTop_solarTide
zTop_gravAcc = zTop_solarTide - 0.50e5
gravAccText = ax1.text( xTop_gravAcc, zTop_gravAcc, '', fontsize=12 )
xGrav = gravAcc_x[data_indices]
yGrav = gravAcc_y[data_indices]
zGrav = gravAcc_z[data_indices]
gravAccData = np.sqrt( xGrav**2 + yGrav**2 + zGrav**2 )
gravAccData = gravAccData.tolist( )
gravAcc_SciNot = []
for index in range( 0, len( gravAccData ) ):
    gravAcc_SciNot.append( "{:.2e}".format( gravAccData[index] ) )

## define the scatter point and the coordinates for which it has to be animated
point, = ax1.plot( [], [], marker='o' )
z = inertial_position_z[ data_indices ].tolist( )
x = inertial_position_x[ data_indices ].tolist( )
coordinate = [ x, z ]
anim = animation.FuncAnimation( fig, animate, init_func=init,
                                frames=len( timeData ), interval=100, blit=False )

ax1.grid(True)
ax1.set_xlabel( 'x [m]' )
ax1.set_ylabel( 'z [m]' )
ax1.set_xlim( min( inertial_position_x ), max( inertial_position_x ) )
ax1.set_ylim( min( inertial_position_z ), max( inertial_position_z ) )
ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
ax1.set_aspect( 1 )

if database:
    database.close( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## show the plot
kwargDict = dict( kw_args = { 'bbox_inches' : 'tight' } )
anim.save( savePath, savefig_kwargs=kwargDict, bitrate=100 )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
