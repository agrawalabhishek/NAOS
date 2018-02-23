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
savePath1 = "../data/guarantee_escape_speed/"
savePath2 = "longest_edge/"
savePath3 = "test.mp4"
savePath = savePath1 + savePath2 + savePath3

try:
    database = sqlite3.connect("../data/guarantee_escape_speed/"
                               + "longest_edge/"
                               + "longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data now...\n"

data = pd.read_sql( "SELECT     ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth ),                                        \
                                ROUND( launch_declination ),                                    \
                                position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                inertial_position_x,                                            \
                                inertial_position_y,                                            \
                                inertial_position_z,                                            \
                                inertial_velocity_x,                                            \
                                inertial_velocity_y,                                            \
                                inertial_velocity_z,                                            \
                                time                                                            \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( launch_azimuth ) = 270.0                                 \
                     AND        ROUND( launch_declination ) = 15.0                              \
                     AND        ROUND( initial_velocity_magnitude ) = 6.0;",                    \
                     database )

if database:
    database.close( )

data.columns = [ 'initial_velocity_magnitude',                                                  \
                 'launch_azimuth',                                                              \
                 'launch_declination',                                                          \
                 'position_x',                                                                  \
                 'position_y',                                                                  \
                 'position_z',                                                                  \
                 'inertial_position_x',                                                         \
                 'inertial_position_y',                                                         \
                 'inertial_position_z',                                                         \
                 'inertial_velocity_x',                                                         \
                 'inertial_velocity_y',                                                         \
                 'inertial_velocity_z',                                                         \
                 'time' ]

initialVelocityMagnitude    = data[ 'initial_velocity_magnitude' ]
launchAzimuth               = data[ 'launch_azimuth' ]
launchDeclination           = data[ 'launch_declination' ]
position_x                  = data[ 'position_x' ]
position_y                  = data[ 'position_y' ]
position_z                  = data[ 'position_z' ]
inertial_position_x         = data[ 'inertial_position_x' ]
inertial_position_y         = data[ 'inertial_position_y' ]
inertial_position_z         = data[ 'inertial_position_z' ]
inertial_velocity_x         = data[ 'inertial_velocity_x' ]
inertial_velocity_y         = data[ 'inertial_velocity_y' ]
inertial_velocity_z         = data[ 'inertial_velocity_z' ]
t                           = data[ 'time' ]

print "Processing data now...\n"

plt.suptitle( '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Launch declination = ' + str( launchDeclination[ 0 ] ) + ' [deg]', fontsize=12 )

# Draw ellipse through the patches package
acwRotationAngle = 0.0
ellipse = mpatches.Ellipse( ( 0.0, 0.0 ), 2.0*alpha, 2.0*beta, acwRotationAngle,
                            alpha=1.0, fill=False,
                            edgecolor='red' )
ax1.add_patch( ellipse )

# Draw the entire trajectory path
ax1.plot( inertial_position_x, inertial_position_y, linewidth=1, color=colors.cnames['purple'] )

def init( ):
    ax1.plot( inertial_position_x[ 0 ], inertial_position_y[ 0 ] )
    return ellipse

def animate( i ):
    timeText.set_text( time_template % ( timeData[ i ] ) )
    velocityText.set_text( velocity_template % ( velocityData[ i ] ) )
    point.set_data( coordinate[ 0 ][ i ], coordinate[ 1 ][ i ] )
    ellipse.angle = RotationAngles[ i ]
    return ellipse, point, timeText, velocityText

### main segment of the code
## data indices
steps = 100
data_indices = range( 0, len( inertial_position_x ), steps )

## define the rotation angles for the ellipse
RotationAngles = ( Wz * t[ data_indices ] ) * 180.0 / np.pi
RotationAngles = RotationAngles.tolist( )

## define the text showing change in time value
time_template = 'Time = %.2f [days]'
# define text placement coordinates
# xTop_text = max( inertial_position_x ) + 0.25e5
# yTop_text = max( inertial_position_y )
xTop_text = 0.0
yTop_text = -0.4e6
timeText = ax1.text( xTop_text, yTop_text, '', fontsize=12 )
timeData = t[ data_indices ] / ( 24.0 * 60.0 * 60.0 )
timeData = timeData.tolist( )

## define the text handle showing velocity magnitude of the regolith
velocity_template = 'Velocity = %0.3f [m/s]'
# define velocity placement coordinates
xTop_velocity = xTop_text
yTop_velocity = yTop_text - 0.50e5
velocityText = ax1.text( xTop_velocity, yTop_velocity, '', fontsize=12 )
Vx = inertial_velocity_x[data_indices]
Vy = inertial_velocity_y[data_indices]
Vz = inertial_velocity_z[data_indices]
velocityData = np.sqrt( Vx**2 + Vy**2 + Vz**2 )
velocityData = velocityData.tolist( )

## define the scatter point and the coordinates for which it has to be animated
point, = ax1.plot( [], [], marker='o' )
y = inertial_position_y[ data_indices ].tolist( )
x = inertial_position_x[ data_indices ].tolist( )
coordinate = [ x, y ]
anim = animation.FuncAnimation( fig, animate, init_func=init,
                                frames=len( RotationAngles ), interval=100, blit=False )

ax1.grid(True)
ax1.set_xlabel( 'x [m]' )
ax1.set_ylabel( 'y [m]' )
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
