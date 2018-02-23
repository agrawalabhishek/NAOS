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

def plotEllipse( semiMajor, semiMinor, angleRange, plotHandle ):
        r = semiMajor * semiMinor / np.sqrt( ( semiMinor * np.cos( angleRange ) )**2
                                           + ( semiMajor * np.sin( angleRange ) )**2 )
        x = r * np.cos( angleRange )
        y = r * np.sin( angleRange )
        plotHandle.plot( x, y )
        plotHandle.grid( )

## Operations
# Connect to SQLite database.
try:
        # database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

# lowerTimeBound = 48.0 * 60.0 * 60.0
# upperTimeBound = 72.0 * 60.0 * 60.0

azimuth_range = tuple( ( np.arange( 0.0, 360.0, 20.0 ) ).tolist( ) )

data1 = pd.read_sql( "SELECT    position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                ROUND( initial_velocity_magnitude ),                        \
                                time,                                                       \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ROUND( initial_velocity_magnitude ) = 10                    \
                     AND        ROUND( launch_azimuth ) IN " + str(azimuth_range) + ";",    \
                     # AND        time >= " + str( lowerTimeBound ) + "                       \
                     # AND        time < " + str( upperTimeBound ) + ";",                     \
                     database )

if database:
    database.close( )

data1.columns = [ 'x',                                                                      \
                  'y',                                                                      \
                  'z',                                                                      \
                  'init_vel_mag',                                                           \
                  'time',                                                                   \
                  'launch_azimuth' ]

x                   = data1[ 'x' ]
y                   = data1[ 'y' ]
z                   = data1[ 'z' ]
initial_velocity    = data1[ 'init_vel_mag' ]
t                   = data1[ 'time' ]
azimuth             = data1[ 'launch_azimuth' ]

# start position for all cases is the same
xPositionStart = x[ 0 ]
yPositionStart = y[ 0 ]
zPositionStart = z[ 0 ]

theta = np.linspace(0, 2 * np.pi, 300)

print "processing data now..."

## for low altitude orbit
fig = plt.figure( )
gs = gridspec.GridSpec( 3, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

plotEllipse( alpha, beta, theta, ax1 )
plotEllipse( beta, gamma, theta, ax2 )
plotEllipse( alpha, gamma, theta, ax3 )

rangeOfParticles = np.sqrt( x**2 + y**2 + z**2 )

LAO = np.where( (rangeOfParticles <= 2.0*alpha) & (rangeOfParticles > 1.0*alpha) )
LAO = LAO[ 0 ]

xPlot = x[ LAO ]
yPlot = y[ LAO ]
zPlot = z[ LAO ]
tPlot = t[ LAO ]
azimuthPlot = azimuth[ LAO ]
unique_azimuths = np.unique( azimuthPlot )
colors = plt.cm.Vega20( np.linspace( 0, 1, len( unique_azimuths ) ) )

for index in range( 0, len( unique_azimuths ) ):
    currentAzimuth = unique_azimuths[ index ]
    dataIndices = np.where( azimuthPlot == currentAzimuth )
    dataIndices = dataIndices[ 0 ]
    xPlotNew = xPlot[ dataIndices ]
    yPlotNew = yPlot[ dataIndices ]
    zPlotNew = zPlot[ dataIndices ]
    ax1.scatter( xPlotNew, yPlotNew, s=2, c=colors[index], edgecolors='face', alpha=1.0, \
                 label='Launch Azimuth = ' + str(currentAzimuth) + '[deg]' )
    ax2.scatter( yPlotNew, zPlotNew, s=2, c=colors[index], edgecolors='face', alpha=1.0 )
    ax3.scatter( xPlotNew, zPlotNew, s=2, c=colors[index], edgecolors='face', alpha=1.0 )

ax1.legend( markerscale=7 ).draggable( )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.grid( True )

ax2.set_xlabel('y [m]')
ax2.set_ylabel('z [m]')
ax2.grid( True )

ax3.set_xlabel('x [m]')
ax3.set_ylabel('z [m]')
ax3.grid( True )


# xPlot = ( x[ LAO ] ).tolist( )
# yPlot = ( y[ LAO ] ).tolist( )
# zPlot = ( z[ LAO ] ).tolist( )
# tPlot = ( t[ LAO ] ).tolist( )


# previousTime = 0.0
# for index in range( 0, len( tPlot ), 10 ):
#     # currentTime = tPlot[ index ]
#     # timeDiff = np.abs( currentTime - previousTime )
#     # if timeDiff > 50:
#     ax1.scatter( xPlot[ index ], yPlot[ index ], s=3, c='red', edgecolors='face', alpha=1.0 )
#     ax2.scatter( yPlot[ index ], zPlot[ index ], s=3, c='red', edgecolors='face', alpha=1.0 )
#     ax3.scatter( xPlot[ index ], zPlot[ index ], s=3, c='red', edgecolors='face', alpha=1.0 )
    # else:
    #     ax1.scatter( xPlot[ index ], yPlot[ index ], s=3, c='green', edgecolors='face', alpha=1.0 )
    #     ax2.scatter( yPlot[ index ], zPlot[ index ], s=3, c='green', edgecolors='face', alpha=1.0 )
    #     ax3.scatter( xPlot[ index ], zPlot[ index ], s=3, c='green', edgecolors='face', alpha=1.0 )
    # previousTime = currentTime


## for medium altitude orbit
# MAO = np.where( (rangeOfParticles > 2.0*alpha) & (rangeOfParticles <= 3.0*alpha) )
# MAO = MAO[ 0 ]
# xPlot = x[ MAO ]
# yPlot = y[ MAO ]
# zPlot = z[ MAO ]

# ax1.scatter( xPlot, yPlot, s=1, c='green', edgecolors='face', alpha=0.1, label='MAO' )
# ax2.scatter( yPlot, zPlot, s=1, c='green', edgecolors='face', alpha=0.1, label='MAO' )
# ax3.scatter( xPlot, zPlot, s=1, c='green', edgecolors='face', alpha=0.1, label='MAO' )

# ax1.legend( markerscale=10 ).draggable( )
# ax1.set_xlabel('x [m]')
# ax1.set_ylabel('y [m]')
# ax1.grid( True )
# # ax1.axis('equal')

# ax2.legend( markerscale=10 ).draggable( )
# ax2.set_xlabel('y [m]')
# ax2.set_ylabel('z [m]')
# ax2.grid( True )
# # ax2.axis('equal')

# ax3.legend( markerscale=10 ).draggable( )
# ax3.set_xlabel('x [m]')
# ax3.set_ylabel('z [m]')
# ax3.grid( True )
# # ax3.axis('equal')

# plt.suptitle( 'Regolith trajectory points at different altitudes \n'
#               + 'LAO - 1.0 Alpha to 2.0 Alpha' + ' MAO - 2.0 Alpha to 3.0 Alpha \n'
#               + '$V_{launch}$ = ' + str( initial_velocity[ 0 ] ) + ' [m/s] \n'
#               + 'Time: ' + str( lowerTimeBound / (60.0*60.0) ) + ' [hrs] - ' + str( upperTimeBound / (60.0*60.0) ) + ' [hrs]' )

# ## for high altitude orbit
# fig = plt.figure( )
# gs = gridspec.GridSpec( 3, 1 )
# ax1 = plt.subplot( gs[ 0 ] )
# ax2 = plt.subplot( gs[ 1 ] )
# ax3 = plt.subplot( gs[ 2 ] )

# plotEllipse( alpha, beta, theta, ax1 )
# plotEllipse( beta, gamma, theta, ax2 )
# plotEllipse( alpha, gamma, theta, ax3 )

# HAO = np.where( (rangeOfParticles > 3.0*alpha) & (rangeOfParticles <= 5.0*alpha) )
# HAO = HAO[ 0 ]
# xPlot = x[ HAO ]
# yPlot = y[ HAO ]
# zPlot = z[ HAO ]

# ax1.scatter( xPlot, yPlot, s=1, c='orange', edgecolors='face', alpha=0.1, label='HAO' )
# ax2.scatter( yPlot, zPlot, s=1, c='orange', edgecolors='face', alpha=0.1, label='HAO' )
# ax3.scatter( xPlot, zPlot, s=1, c='orange', edgecolors='face', alpha=0.1, label='HAO' )

# ax1.legend( markerscale=10 ).draggable( )
# ax1.set_xlabel('x [m]')
# ax1.set_ylabel('y [m]')
# ax1.grid( True )
# # ax1.axis('equal')

# ax2.legend( markerscale=10 ).draggable( )
# ax2.set_xlabel('y [m]')
# ax2.set_ylabel('z [m]')
# ax2.grid( True )
# # ax2.axis('equal')

# ax3.legend( markerscale=10 ).draggable( )
# ax3.set_xlabel('x [m]')
# ax3.set_ylabel('z [m]')
# ax3.grid( True )
# # ax3.axis('equal')

# plt.suptitle( 'Regolith trajectory points at different altitudes \n'
#               + 'HAO - 3.0 Alpha to 5.0 Alpha \n'
#               + '$V_{launch}$ = ' + str( initial_velocity[ 0 ] ) + ' [m/s] \n'
#               + 'Time: ' + str( lowerTimeBound / (60.0*60.0) ) + ' [hrs] - ' + str( upperTimeBound / (60.0*60.0) ) + ' [hrs]' )

# ## for ultra-high altitude orbit
# fig = plt.figure( )
# gs = gridspec.GridSpec( 3, 1 )
# ax1 = plt.subplot( gs[ 0 ] )
# ax2 = plt.subplot( gs[ 1 ] )
# ax3 = plt.subplot( gs[ 2 ] )

# plotEllipse( alpha, beta, theta, ax1 )
# plotEllipse( beta, gamma, theta, ax2 )
# plotEllipse( alpha, gamma, theta, ax3 )

# UHAO = np.where( (rangeOfParticles > 5.0*alpha) & (rangeOfParticles <= 7.0*alpha) )
# UHAO = UHAO[ 0 ]
# xPlot = x[ UHAO ]
# yPlot = y[ UHAO ]
# zPlot = z[ UHAO ]

# ax1.scatter( xPlot, yPlot, s=1, c=colors.cnames['purple'], edgecolors='face', alpha=0.1, label='UHAO' )
# ax2.scatter( yPlot, zPlot, s=1, c=colors.cnames['purple'], edgecolors='face', alpha=0.1, label='UHAO' )
# ax3.scatter( xPlot, zPlot, s=1, c=colors.cnames['purple'], edgecolors='face', alpha=0.1, label='UHAO' )

# plt.suptitle( 'Regolith trajectory points at different altitudes \n'
#               + 'UHAO - 5.0 Alpha to 7.0 Alpha \n'
#               + '$V_{launch}$ = ' + str( initial_velocity[ 0 ] ) + ' [m/s] \n'
#               + 'Time: ' + str( lowerTimeBound / (60.0*60.0) ) + ' [hrs] - ' + str( upperTimeBound / (60.0*60.0) ) + ' [hrs]' )

# ax1.legend( markerscale=10 ).draggable( )
# ax1.set_xlabel('x [m]')
# ax1.set_ylabel('y [m]')
# ax1.grid( True )
# # ax1.axis('equal')

# ax2.legend( markerscale=10 ).draggable( )
# ax2.set_xlabel('y [m]')
# ax2.set_ylabel('z [m]')
# ax2.grid( True )
# # ax2.axis('equal')

# ax3.legend( markerscale=10 ).draggable( )
# ax3.set_xlabel('x [m]')
# ax3.set_ylabel('z [m]')
# ax3.grid( True )
# # ax3.axis('equal')

# ## for extremely-high altitude orbit
# fig = plt.figure( )
# gs = gridspec.GridSpec( 3, 1 )
# ax1 = plt.subplot( gs[ 0 ] )
# ax2 = plt.subplot( gs[ 1 ] )
# ax3 = plt.subplot( gs[ 2 ] )

# plotEllipse( alpha, beta, theta, ax1 )
# plotEllipse( beta, gamma, theta, ax2 )
# plotEllipse( alpha, gamma, theta, ax3 )

# EHAO = np.where( (rangeOfParticles > 7.0*alpha) )
# EHAO = EHAO[ 0 ]
# xPlot = x[ EHAO ]
# yPlot = y[ EHAO ]
# zPlot = z[ EHAO ]

# ax1.scatter( xPlot, yPlot, s=1, c=colors.cnames['purple'], edgecolors='face', alpha=0.1, label='EHAO' )
# ax2.scatter( yPlot, zPlot, s=1, c=colors.cnames['purple'], edgecolors='face', alpha=0.1, label='EHAO' )
# ax3.scatter( xPlot, zPlot, s=1, c=colors.cnames['purple'], edgecolors='face', alpha=0.1, label='EHAO' )

# plt.suptitle( 'Regolith trajectory points at different altitudes \n'
#               + 'EHAO - 7.0 Alpha and above \n'
#               + '$V_{launch}$ = ' + str( initial_velocity[ 0 ] ) + ' [m/s] \n'
#               + 'Time: ' + str( lowerTimeBound / (60.0*60.0) ) + ' [hrs] - ' + str( upperTimeBound / (60.0*60.0) ) + ' [hrs]' )

# ax1.legend( markerscale=10 ).draggable( )
# ax1.set_xlabel('x [m]')
# ax1.set_ylabel('y [m]')
# ax1.grid( True )
# # ax1.axis('equal')

# ax2.legend( markerscale=10 ).draggable( )
# ax2.set_xlabel('y [m]')
# ax2.set_ylabel('z [m]')
# ax2.grid( True )
# # ax2.axis('equal')

# ax3.legend( markerscale=10 ).draggable( )
# ax3.set_xlabel('x [m]')
# ax3.set_ylabel('z [m]')
# ax3.grid( True )
# ax3.axis('equal')

if database:
    database.close( )

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
