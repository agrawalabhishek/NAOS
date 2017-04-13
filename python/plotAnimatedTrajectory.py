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
                 'time' ]

velocityMagnitude   = data[ 'initial_velocity_magnitude' ]
launchAzimuth       = data[ 'launch_azimuth' ]
solarPhase          = data[ 'initial_solar_phase_angle' ]
changingSolarPhase  = data[ 'solar_phase_angle' ]
position_x          = data[ 'position_x' ]
position_y          = data[ 'position_y' ]
position_z          = data[ 'position_z' ]
inertial_position_x = data[ 'inertial_position_x' ]
inertial_position_y = data[ 'inertial_position_y' ]
inertial_position_z = data[ 'inertial_position_z' ]
t                   = data[ 'time' ]

print "Processing data now...\n"

print "Solar phase at start of the segment = " + str( changingSolarPhase[ 0 ] ) + " [deg]"

# Set up the figure
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

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
    point.set_data( coordinate[ 0 ][ i ], coordinate[ 1 ][ i ] )
    ellipse.angle = RotationAngles[ i ]
    return ellipse, point

## main segment of the code
# data indices
steps = 10
data_indices = range( 0, len( inertial_position_x ), steps )

# define the rotation angles for the ellipse
RotationAngles = ( Wz * t[ data_indices ] ) * 180.0 / np.pi
RotationAngles = RotationAngles.tolist( )

# define the scatter point and the coordinates for which it has to be animated
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
anim.save(savePath)
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
