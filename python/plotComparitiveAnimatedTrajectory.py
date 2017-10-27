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
fig = plt.figure( figsize=(20, 20) )
gs = gridspec.GridSpec( 1, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )

# Connect to SQLite database.
savePath1 = "../data/regolith_launched_from_longest_edge/"
savePath2 = "multiple_launch_velocity_with_perturbations/"
savePath3 = "simulation_time_9_months/"
savePath4 = "3.2Density_1cmSize/"
savePath5 = "test.mp4"
savePath = savePath1 + savePath2 + savePath3 + savePath4 + savePath5

## animation for the case without solar perturbations
try:
       database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                                   + "multiple_launch_velocity/"
                                   + "simulation_time_9_months/"
                                   + "longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data for the non-perturbing case now...\n"

data = pd.read_sql( "SELECT     position_x,                                         \
                                position_y,                                         \
                                position_z,                                         \
                                ROUND( initial_velocity_magnitude ),                \
                                inertial_position_x,                                \
                                inertial_position_y,                                \
                                inertial_position_z,                                \
                                inertial_velocity_x,                                \
                                inertial_velocity_y,                                \
                                inertial_velocity_z,                                \
                                ROUND( launch_azimuth ),                            \
                                time                                                \
                     FROM       regolith_trajectory_results                         \
                     WHERE      ROUND( launch_azimuth ) = 165.0                     \
                     AND        ROUND( initial_velocity_magnitude ) = 8.0;",        \
                     database )

data.columns = [ 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'velocity_magnitude',                                  \
                 'inertial_position_x',                                 \
                 'inertial_position_y',                                 \
                 'inertial_position_z',                                 \
                 'inertial_velocity_x',                                 \
                 'inertial_velocity_y',                                 \
                 'inertial_velocity_z',                                 \
                 'launch_azimuth',                                      \
                 'time' ]

x                   = data[ 'x' ]
y                   = data[ 'y' ]
z                   = data[ 'z' ]
velocityMagnitude   = data[ 'velocity_magnitude' ]
inertial_x          = data[ 'inertial_position_x' ]
inertial_y          = data[ 'inertial_position_y' ]
inertial_z          = data[ 'inertial_position_z' ]
inertial_velocity_x = data[ 'inertial_velocity_x' ]
inertial_velocity_y = data[ 'inertial_velocity_y' ]
inertial_velocity_z = data[ 'inertial_velocity_z' ]
launchAzimuth       = data[ 'launch_azimuth' ]
t                   = data[ 'time' ]

if database:
    database.close( )

# get the end point for a data array
endIndex = np.size( x )

# Draw ellipse through the patches package
acwRotationAngle = 0.0
ellipse_noSolarPerturbations = mpatches.Ellipse( ( 0.0, 0.0 ), 2.0*alpha, 2.0*beta, acwRotationAngle,
                                                 alpha=1.0, fill=False,
                                                 edgecolor='red' )
ax1.add_patch( ellipse_noSolarPerturbations )

# Draw the entire trajectory path
ax1.plot( inertial_x, inertial_y, linewidth=1, color=colors.cnames['purple'] )

## data indices
steps = 10
data_indices = range( 0, len( inertial_x ), steps )

## define the rotation angles for the ellipse
RotationAngles_noSolarPerturbations = ( Wz * t[ data_indices ] ) * 180.0 / np.pi
RotationAngles_noSolarPerturbations = RotationAngles_noSolarPerturbations.tolist( )

## define the text showing change in time value
time_template = 'Time = %.2f [hrs]'
# define text placement coordinates
# xTop_text = max( inertial_position_x ) + 0.25e5
# yTop_text = max( inertial_position_y )
xTop_text = 0.5e5
yTop_text = 0.6e5
nextLineAt = 0.1e5
timeText_noSolarPerturbations = ax1.text( xTop_text, yTop_text, '', fontsize=12 )
timeData = t[ data_indices ] / ( 60.0 * 60.0 )
timeData_noSolarPerturbations = timeData.tolist( )

## define the text handle showing velocity magnitude of the regolith
velocity_template = 'Velocity = %0.3f [m/s]'
# define velocity placement coordinates
xTop_velocity = xTop_text
yTop_velocity = yTop_text - nextLineAt
velocityText_noSolarPerturbations = ax1.text( xTop_velocity, yTop_velocity, '', fontsize=12 )
Vx = inertial_velocity_x[data_indices]
Vy = inertial_velocity_y[data_indices]
Vz = inertial_velocity_z[data_indices]
velocityData = np.sqrt( Vx**2 + Vy**2 + Vz**2 )
velocityData_noSolarPerturbations = velocityData.tolist( )

## define the text handle showing range magnitude of the regolith
range_template = 'Range = %0.3f [m]'
# define velocity placement coordinates
xTop_range = xTop_velocity
yTop_range = yTop_velocity - nextLineAt
rangeText_noSolarPerturbations = ax1.text( xTop_range, yTop_range, '', fontsize=12 )
rangeData = np.sqrt( inertial_x[data_indices]**2 + inertial_y[data_indices]**2 + inertial_z[data_indices]**2 )
rangeData_noSolarPerturbations = rangeData.tolist( )

## define the scatter point and the coordinates for which it has to be animated
point_noSolarPerturbations, = ax1.plot( [], [], marker='o' )
y_noSolarPerturbations = inertial_y[ data_indices ].tolist( )
x_noSolarPerturbations = inertial_x[ data_indices ].tolist( )
coordinate_noSolarPerturbations = [ x_noSolarPerturbations, y_noSolarPerturbations ]

if database:
    database.close( )

## animation for the case with perturbations from the Sun
try:
    database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                               + "multiple_launch_velocity_with_perturbations/"
                               + "simulation_time_9_months/"
                               + "3.2Density_1cmSize/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data for the perturbing case now...\n"

lowerTimeDays = 0.0
lowerTime = lowerTimeDays * 24.0 * 60.0 * 60.0

upperTime = t[endIndex-1]
upperTimeDays = upperTime / ( 24.0 * 60.0 * 60.0 )

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
                     WHERE      ROUND( launch_azimuth ) = 165.0                                 \
                     AND        ROUND( initial_velocity_magnitude ) = 8.0                       \
                     AND        ROUND( initial_solar_phase_angle ) = 45.0                       \
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

t_solarPerturbations        = data[ 'time' ]

print "Processing data now...\n"

print "Solar phase at start of the segment = " + str( changingSolarPhase[ 0 ] ) + " [deg]"

plt.suptitle( '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]', fontsize=12 )

asteroidRotationPeriod = ( 2.0 * np.pi / Wz ) / ( 60.0 * 60.0 )
print "Asteroid rotation period = " + str( asteroidRotationPeriod ) + " [hrs]\n"

# Draw ellipse through the patches package
acwRotationAngle = 0.0
ellipse = mpatches.Ellipse( ( 0.0, 0.0 ), 2.0*alpha, 2.0*beta, acwRotationAngle,
                            alpha=1.0, fill=False,
                            edgecolor='red' )
ax2.add_patch( ellipse )

# Draw the entire trajectory path
ax2.plot( inertial_position_x, inertial_position_y, linewidth=1, color=colors.cnames['purple'] )

def init( ):
    ax1.plot( inertial_x[ 0 ], inertial_y[ 0 ] )
    ax2.plot( inertial_position_x[ 0 ], inertial_position_y[ 0 ] )
    return ellipse, ellipse_noSolarPerturbations

def animate( i ):
    # non perturbing data
    timeText_noSolarPerturbations.set_text( time_template % ( timeData_noSolarPerturbations[ i ] ) )
    velocityText_noSolarPerturbations.set_text( velocity_template % ( velocityData_noSolarPerturbations[ i ] ) )
    rangeText_noSolarPerturbations.set_text( range_template % ( rangeData_noSolarPerturbations[ i ] ) )
    point_noSolarPerturbations.set_data( coordinate_noSolarPerturbations[ 0 ][ i ], coordinate_noSolarPerturbations[ 1 ][ i ] )
    ellipse_noSolarPerturbations.angle = RotationAngles_noSolarPerturbations[ i ]

    # perturbing data
    timeText.set_text( time_template % ( timeData[ i ] ) )
    velocityText.set_text( velocity_template % ( velocityData[ i ] ) )
    rangeText.set_text( range_template % ( rangeData[ i ] ) )
    srpText.set_text( srp_template % ( srpData_SciNot[ i ] ) )
    solarTideText.set_text( solarTide_template % ( solarTideData_SciNot[ i ] ) )
    gravAccText.set_text( gravAcc_template % ( gravAcc_SciNot[ i ] ) )
    point.set_data( coordinate[ 0 ][ i ], coordinate[ 1 ][ i ] )
    ellipse.angle = RotationAngles[ i ]

    return timeText_noSolarPerturbations, velocityText_noSolarPerturbations,                    \
           rangeText_noSolarPerturbations,                                                      \
           point_noSolarPerturbations, ellipse_noSolarPerturbations,                            \
           ellipse, point, timeText, velocityText, rangeText, srpText, solarTideText, gravAccText

##perturbing case data processing
## data indices
# steps = 10
# data_indices = range( 0, len( inertial_position_x ), steps )

timeToExtract =  t[data_indices].tolist( )
data_indices = []
for index in range( 0, len( timeToExtract ) ):
    currentIndex = np.where( t_solarPerturbations == timeToExtract[ index ] )
    currentIndex = currentIndex[ 0 ]
    if currentIndex.size != 0:
        data_indices.append( int( currentIndex ) )

timeValuesExtracted = t_solarPerturbations[ data_indices ]

## define the rotation angles for the ellipse
RotationAngles = ( Wz * t_solarPerturbations[ data_indices ] ) * 180.0 / np.pi
RotationAngles = RotationAngles.tolist( )

## define the text showing change in time value
time_template = 'Time = %.2f [hrs]'
# define text placement coordinates
# xTop_text = max( inertial_position_x ) + 0.25e5
# yTop_text = max( inertial_position_y )
xTop_text = 0.4e5
yTop_text = 6.0e4
timeText = ax2.text( xTop_text, yTop_text, '', fontsize=12 )
timeData = t_solarPerturbations[ data_indices ] / ( 60.0 * 60.0 )
timeData = timeData.tolist( )

nextLineAt = 1.0e4
## define the text handle showing velocity magnitude of the regolith
velocity_template = 'Velocity = %0.3f [m/s]'
# define velocity placement coordinates
xTop_velocity = xTop_text
yTop_velocity = yTop_text - nextLineAt
velocityText = ax2.text( xTop_velocity, yTop_velocity, '', fontsize=12 )
Vx = inertial_velocity_x[data_indices]
Vy = inertial_velocity_y[data_indices]
Vz = inertial_velocity_z[data_indices]
velocityData = np.sqrt( Vx**2 + Vy**2 + Vz**2 )
velocityData = velocityData.tolist( )

## define the text handle showing range magnitude of the regolith
range_template = 'Range = %0.3f [m]'
# define velocity placement coordinates
xTop_range = xTop_velocity
yTop_range = yTop_velocity - nextLineAt
rangeText = ax2.text( xTop_range, yTop_range, '', fontsize=12 )
rangeData = np.sqrt( inertial_position_x[data_indices]**2 + inertial_position_y[data_indices]**2 + inertial_position_z[data_indices]**2 )
rangeData = rangeData.tolist( )

## define the text handle showing various accelerations acting on the regolith
# SRP
srp_template = 'SRP = %s $[m/s^2]$'
xTop_srp = xTop_range
yTop_srp = yTop_range - nextLineAt
srpText = ax2.text( xTop_srp, yTop_srp, '', fontsize=12 )
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
yTop_solarTide = yTop_srp - nextLineAt
solarTideText = ax2.text( xTop_solarTide, yTop_solarTide, '', fontsize=12 )
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
yTop_gravAcc = yTop_solarTide - nextLineAt
gravAccText = ax2.text( xTop_gravAcc, yTop_gravAcc, '', fontsize=12 )
xGrav = gravAcc_x[data_indices]
yGrav = gravAcc_y[data_indices]
zGrav = gravAcc_z[data_indices]
gravAccData = np.sqrt( xGrav**2 + yGrav**2 + zGrav**2 )
gravAccData = gravAccData.tolist( )
gravAcc_SciNot = []
for index in range( 0, len( gravAccData ) ):
    gravAcc_SciNot.append( "{:.2e}".format( gravAccData[index] ) )

## define the scatter point and the coordinates for which it has to be animated
point, = ax2.plot( [], [], marker='o' )
y = inertial_position_y[ data_indices ].tolist( )
x = inertial_position_x[ data_indices ].tolist( )
coordinate = [ x, y ]
anim = animation.FuncAnimation( fig, animate, init_func=init,
                                frames=len( RotationAngles ), interval=100, blit=False )

ax1.grid(True)
ax1.set_xlabel( 'x [m]' )
ax1.set_ylabel( 'y [m]' )
# ax1.set_xlim( min( inertial_x ), 1.0e5 )
# ax1.set_ylim( min( inertial_position_y ), 0.4e5 )
ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
ax1.set_aspect( 1 )
ax1.set_title( 'No Solar Perturbations' )

ax2.grid(True)
ax2.set_xlabel( 'x [m]' )
ax2.set_ylabel( 'y [m]' )
# ax2.set_xlim( min( inertial_x ), 1.0e5 )
# ax2.set_ylim( min( inertial_position_y ), 0.4e5 )
ax2.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
ax2.set_aspect( 1 )
ax2.set_title( 'SRP + STBE' )

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
