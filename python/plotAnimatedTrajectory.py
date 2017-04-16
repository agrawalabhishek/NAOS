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

upperTimeDays = 10.0
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
                     WHERE      ROUND( launch_azimuth ) = 185.0                                 \
                     AND        ROUND( initial_velocity_magnitude ) = 5.0                       \
                     AND        ROUND( initial_solar_phase_angle ) = 135.0                      \
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

# Set up the figure
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
plt.suptitle( '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]', fontsize=10 )

# Draw ellipse through the patches package
acwRotationAngle = 0.0
ellipse = mpatches.Ellipse( ( 0.0, 0.0 ), 2.0*alpha, 2.0*beta, acwRotationAngle,
                            alpha=1.0, fill=False,
                            edgecolor='red' )
ax1.add_patch( ellipse )

# Draw the entire trajectory path
ax1.plot( inertial_position_x, inertial_position_y, linewidth=1, color=colors.cnames['purple'] )

# mark starting sun position as a quiver arrow
# startSunAngle = changingSolarPhase[ 0 ] * np.pi / 180.0
# sunStartPos_x = 1.0e10 * np.cos( startSunAngle )
# sunStartPos_y = 1.0e10 * np.sin( startSunAngle )
# ax1.quiver( 0.0, 0.0, sunStartPos_x, sunStartPos_y, pivot='tail', color='magenta' )

# mark end position of sun as a quiver arrow
# endSunAngle = changingSolarPhase[ len( changingSolarPhase ) - 1 ] * np.pi / 180.0
# sunEndPos_x = 1.0e10 * np.cos( endSunAngle )
# sunEndPos_y = 1.0e10 * np.sin( endSunAngle )
# ax1.quiver( 0.0, 0.0, sunEndPos_x, sunEndPos_y, pivot='tail', color='black' )

def init( ):
    ax1.plot( inertial_position_x[ 0 ], inertial_position_y[ 0 ] )
    return ellipse

def animate( i ):
    timeText.set_text( time_template % ( timeData[ i ] ) )
    velocityText.set_text( velocity_template % ( velocityData[ i ] ) )
    srpText.set_text( srp_template % ( srpData_SciNot[ i ] ) )
    solarTideText.set_text( solarTide_template % ( solarTideData_SciNot[ i ] ) )
    gravAccText.set_text( gravAcc_template % ( gravAcc_SciNot[ i ] ) )
    point.set_data( coordinate[ 0 ][ i ], coordinate[ 1 ][ i ] )
    ellipse.angle = RotationAngles[ i ]
    return ellipse, point, timeText, velocityText, srpText, solarTideText, gravAccText

### main segment of the code
## data indices
steps = 50
data_indices = range( 0, len( inertial_position_x ), steps )

## define the rotation angles for the ellipse
RotationAngles = ( Wz * t[ data_indices ] ) * 180.0 / np.pi
RotationAngles = RotationAngles.tolist( )

## define the text showing change in time value
time_template = 'Time = %.2f [days]'
# define text placement coordinates
xTop_text = max( inertial_position_x ) + 0.25e5
yTop_text = max( inertial_position_y )
timeText = ax1.text( xTop_text, yTop_text, '', fontsize=8 )
timeData = t[ data_indices ] / ( 24.0 * 60.0 * 60.0 )
timeData = timeData.tolist( )

## define the text handle showing velocity magnitude of the regolith
velocity_template = 'Velocity = %0.2f [m/s]'
# define velocity placement coordinates
xTop_velocity = xTop_text
yTop_velocity = yTop_text - 0.25e5
velocityText = ax1.text( xTop_velocity, yTop_velocity, '', fontsize=8 )
velocityData = np.sqrt( inertial_velocity_x**2 + inertial_velocity_y**2 + inertial_velocity_z**2 )
velocityData = velocityData.tolist( )

## define the text handle showing various accelerations acting on the regolith
# SRP
srp_template = 'SRP = %s $[m^2/s^2]$'
xTop_srp = xTop_velocity
yTop_srp = yTop_velocity - 0.25e5
srpText = ax1.text( xTop_srp, yTop_srp, '', fontsize=8 )
srpData = np.sqrt( srp_x**2 + srp_y**2 + srp_z**2 )
srpData = srpData.tolist( )
srpData_SciNot = []
for index in range( 0, len( srpData ) ):
    srpData_SciNot.append( "{:.2e}".format( srpData[index] ) )

# Solar tide
solarTide_template = 'Solar tide = %s $[m^2/s^2]$'
xTop_solarTide = xTop_srp
yTop_solarTide = yTop_srp - 0.25e5
solarTideText = ax1.text( xTop_solarTide, yTop_solarTide, '', fontsize=8 )
solarTideData = np.sqrt( solarTide_x**2 + solarTide_y**2 + solarTide_z**2 )
solarTideData = solarTideData.tolist( )
solarTideData_SciNot = []
for index in range( 0, len( solarTideData ) ):
    solarTideData_SciNot.append( "{:.2e}".format( solarTideData[index] ) )

# Gravity
gravAcc_template = 'Grav. = %s $[m^2/s^2]$'
xTop_gravAcc = xTop_solarTide
yTop_gravAcc = yTop_solarTide - 0.25e5
gravAccText = ax1.text( xTop_gravAcc, yTop_gravAcc, '', fontsize=8 )
gravAccData = np.sqrt( gravAcc_x**2 + gravAcc_y**2 + gravAcc_z**2 )
gravAccData = gravAccData.tolist( )
gravAcc_SciNot = []
for index in range( 0, len( gravAccData ) ):
    gravAcc_SciNot.append( "{:.2e}".format( gravAccData[index] ) )

## define the scatter point and the coordinates for which it has to be animated
point, = ax1.plot( [], [], marker='o' )
y = inertial_position_y[ data_indices ].tolist( )
x = inertial_position_x[ data_indices ].tolist( )
coordinate = [ x, y ]
anim = animation.FuncAnimation( fig, animate, init_func=init,
                                frames=len( RotationAngles ), interval=100, blit=False )

ax1.grid(True)
ax1.set_xlim( min( inertial_position_x ), max( inertial_position_x ) )
ax1.set_ylim( min( inertial_position_y ), max( inertial_position_y ) )
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
