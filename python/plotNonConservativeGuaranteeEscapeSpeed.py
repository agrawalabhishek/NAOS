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
Wz = 0.00033118202125129593
mu = 876514

## Set up the figure
fig = plt.figure( )
gs = gridspec.GridSpec( 2, 1, height_ratios = [ 1, 1 ] )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
# ax1 = fig.add_subplot(111)
# ax2 = fig.add_subplot(212)

## Operations
# Connect to SQLite database.
try:
        database = sqlite3.connect("../data/guarantee_escape_speed/longest_edge/longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

data1 = pd.read_sql( "SELECT    velocity_x,                                                 \
                                velocity_y,                                                 \
                                velocity_z,                                                 \
                                inertial_velocity_x,                                        \
                                inertial_velocity_y,                                        \
                                inertial_velocity_z,                                        \
                                time,                                                       \
                                directional_escape_speed,                                   \
                                directional_inertial_escape_speed,                          \
                                ROUND( launch_declination )                                 \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( start_flag = 1 );",                                       \
                     database )

data1.columns = [ 'rotFrame_vx',                                                            \
                  'rotFrame_vy',                                                            \
                  'rotFrame_vz',                                                            \
                  'inertial_vx',                                                            \
                  'inertial_vy',                                                            \
                  'inertial_vz',                                                            \
                  'time',                                                                   \
                  'directional_escape_speed',                                               \
                  'directional_inertial_escape_speed',                                      \
                  'launch_declination' ]

rotFrame_vx                     = data1[ 'rotFrame_vx' ]
rotFrame_vy                     = data1[ 'rotFrame_vy' ]
rotFrame_vz                     = data1[ 'rotFrame_vz' ]
inertial_vx                     = data1[ 'inertial_vx' ]
inertial_vy                     = data1[ 'inertial_vy' ]
inertial_vz                     = data1[ 'inertial_vz' ]
t                               = data1[ 'time' ]
directionalEscapeSpeed          = data1[ 'directional_escape_speed' ]
inertialDirectionalEscapeSpeed  = data1[ 'directional_inertial_escape_speed' ]
declination                     = data1[ 'launch_declination' ]

## get azimuth values for escape and reimpact cases
data2 = pd.read_sql( "SELECT    ROUND( launch_declination )                                 \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 );",                                      \
                     database )

data2.columns = [ 'escape_declination' ]
escape_declination = data2[ 'escape_declination' ]

data3 = pd.read_sql( "SELECT    ROUND( launch_declination )                                 \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( crash_flag = 1 );",                                       \
                     database )

data3.columns = [ 'crash_declination' ]
crash_declination = data3[ 'crash_declination' ]

## extract data from csv file for the new escape speed values
data4 = pd.read_csv( "../data/guarantee_escape_speed/longest_edge/no_solar_perturbations/nonConventionalEscape_qinfinity_0.1alpha.csv" )
data5 = pd.read_csv( "../data/guarantee_escape_speed/longest_edge/no_solar_perturbations/nonConventionalEscape_qinfinity_0.2alpha.csv" )
data6 = pd.read_csv( "../data/guarantee_escape_speed/longest_edge/no_solar_perturbations/nonConventionalEscape_qinfinity_0.3alpha.csv" )

launchDeclination = data4['launch_declination'].values
# print launchDeclination

qInfinity4 = data4['q_infinity'].values
# print qInfinity4
qInfinity5 = data5['q_infinity'].values
qInfinity6 = data6['q_infinity'].values

inertialFrameEscapeSpeedWithPlus4 = data4["inertial_frame_escape_speed_plus"].values
inertialFrameEscapeSpeedWithPlus5 = data5["inertial_frame_escape_speed_plus"].values
inertialFrameEscapeSpeedWithPlus6 = data6["inertial_frame_escape_speed_plus"].values
# print inertialFrameEscapeSpeedWithPlus4

inertialFrameEscapeSpeedWithMinus4 = data4["inertial_frame_escape_speed_minus"].values
inertialFrameEscapeSpeedWithMinus5 = data5["inertial_frame_escape_speed_minus"].values
inertialFrameEscapeSpeedWithMinus6 = data6["inertial_frame_escape_speed_minus"].values

bodyFrameEscapeSpeedWithPlus4 = data4["body_frame_escape_speed_plus"].values
bodyFrameEscapeSpeedWithPlus5 = data5["body_frame_escape_speed_plus"].values
bodyFrameEscapeSpeedWithPlus6 = data6["body_frame_escape_speed_plus"].values
# print bodyFrameEscapeSpeedWithPlus

bodyFrameEscapeSpeedWithMinus4 = data4["body_frame_escape_speed_minus"].values
bodyFrameEscapeSpeedWithMinus5 = data5["body_frame_escape_speed_minus"].values
bodyFrameEscapeSpeedWithMinus6 = data6["body_frame_escape_speed_minus"].values
print bodyFrameEscapeSpeedWithMinus4
print bodyFrameEscapeSpeedWithMinus5
print bodyFrameEscapeSpeedWithMinus6

## plot launch azimuth versus rot frame initial velocity
rotFrameInitialVelocity = np.sqrt( rotFrame_vx**2 + rotFrame_vy**2 + rotFrame_vz**2 )
ax1Handle1 = ax1.plot( declination, rotFrameInitialVelocity, color=colors.cnames['purple'], label='launch speed' )
ax1Handle2 = ax1.plot( declination, directionalEscapeSpeed, color=colors.cnames['black'], label='conservative escape speed' )
ax1.plot( launchDeclination, bodyFrameEscapeSpeedWithPlus4, color=colors.cnames['red'],
    label='$q_\infty$ = ' + str( qInfinity4[ 0 ] / alpha ) + ' alpha' )
ax1.plot( launchDeclination, bodyFrameEscapeSpeedWithPlus5, color=colors.cnames['blue'],
    label='$q_\infty$ = ' + str( qInfinity5[ 0 ] / alpha ) + ' alpha' )
ax1.plot( launchDeclination, bodyFrameEscapeSpeedWithPlus6, color=colors.cnames['green'],
    label='$q_\infty$ = ' + str( qInfinity6[ 0 ] / alpha ) + ' alpha' )

# meanEscapeSpeed = np.mean( [ directionalEscapeSpeed, bodyFrameEscapeSpeedWithPlus4, bodyFrameEscapeSpeedWithPlus5, bodyFrameEscapeSpeedWithPlus6 ], axis=0 )
# ax1.plot( launchDeclination, meanEscapeSpeed, color=colors.cnames['brown'], linestyle='dotted', label='mean escape speed' )

## format axis and title
ax1.set_xlabel('Launch declination [deg]')
ax1.set_ylabel('V [m/s]')
ax1.set_xlim( 10, 80 )
ax1.set_title( 'Guaranteed escape velocity versus launch declination (rotating frame)' )
ax1.legend( ).draggable( )
# ax1.legend( loc='center left', bbox_to_anchor=( 0.995,0.5 ) )
ax1.grid( )

## Plot launch azimuth versus inertial initial velocity
inertialInitialVelocity = np.sqrt( inertial_vx**2 + inertial_vy**2 + inertial_vz**2 )
ax2Handle1, = ax2.plot( declination, inertialInitialVelocity, color=colors.cnames['purple'], label='launch speed' )
ax2Handle2, = ax2.plot( declination, inertialDirectionalEscapeSpeed, color=colors.cnames['black'], label='conservative escape speed' )
# ax2Handle3, = ax2.plot( declination, inertialFrameEscapeSpeedWithPlus, color=colors.cnames['red'], label='non-conservative escape speed' )
ax2Handle4, = ax2.plot( launchDeclination, inertialFrameEscapeSpeedWithPlus4, color=colors.cnames['red'],
    label='$q_\infty$ = ' + str( qInfinity4[ 0 ] / alpha ) + ' alpha' )
ax2Handle5, = ax2.plot( launchDeclination, inertialFrameEscapeSpeedWithPlus5, color=colors.cnames['blue'],
    label='$q_\infty$ = ' + str( qInfinity5[ 0 ] / alpha ) + ' alpha' )
ax2Handle6, = ax2.plot( launchDeclination, inertialFrameEscapeSpeedWithPlus6, color=colors.cnames['green'],
    label='$q_\infty$ = ' + str( qInfinity6[ 0 ] / alpha ) + ' alpha' )


# Find bounds for azimuths that lead to escape
escapeBoundMins = []
escapeBoundMaxs = []

escapeBoundMins.append( min( escape_declination ) )
for index in tqdm( range( len( escape_declination ) ) ):
    if escape_declination[ index ] != max( escape_declination ):
        if escape_declination[ index+1 ] - escape_declination[ index ] != 1:
            escapeBoundMins.append( escape_declination[ index+1 ] )
            escapeBoundMaxs.append( escape_declination[ index ] )
escapeBoundMaxs.append( max( escape_declination ) )

for index in range( len( escapeBoundMins ) ):
    escapeBoundHandle = ax2.axvspan( escapeBoundMins[ index ], escapeBoundMaxs[ index ], color='red', alpha=0.2, label='escape' )
    ax1EscapeHandle = ax1.axvspan( escapeBoundMins[ index ], escapeBoundMaxs[ index ], color='red', alpha=0.2, label='escape' )

# find bounds to azimuth that lead to crash
crashBoundMins = []
crashBoundMaxs = []

crashBoundMins.append( min( crash_declination ) )
for index in tqdm( range( len( crash_declination ) ) ):
    if crash_declination[ index ] != max( crash_declination ):
        if crash_declination[ index+1 ] - crash_declination[ index ] != 1:
            crashBoundMins.append( crash_declination[ index+1 ] )
            crashBoundMaxs.append( crash_declination[ index ] )
crashBoundMaxs.append( max( crash_declination ) )

for index in range( len( crashBoundMins ) ):
    crashBoundHandle = ax2.axvspan( crashBoundMins[ index ], crashBoundMaxs[ index ], color='green', alpha=0.2, label='reimpact' )
    ax1CrashHandle = ax1.axvspan( crashBoundMins[ index ], crashBoundMaxs[ index ], color='green', alpha=0.2, label='reimpact' )

## format axis and title
ax2.set_xlabel('Launch declination [deg]')
ax2.set_ylabel('V [m/s]')
ax2.set_xlim( 10, 80 )
ax2.set_title( 'Guaranteed escape velocity versus launch declination (inertial frame)' )
# ax2.legend( )
# ax2.legend( loc='center left', bbox_to_anchor=( 0.995,0.5 ) )
ax2.grid( )

## Show the plot
plt.legend( handles = [ ax2Handle1, ax2Handle2, ax2Handle4, ax2Handle5, ax2Handle6, escapeBoundHandle, crashBoundHandle ],
            bbox_to_anchor = ( 0.95, 0.5 ),
            loc = 1 )

# plt.legend( handles = [ ax2Handle1, ax2Handle2, escapeBoundHandle ],
#             bbox_to_anchor = ( 0.95, 0.5 ),
#             loc = 1 )

# plt.legend( bbox_to_anchor=(1.0, 1), loc=1, borderaxespad=0. )
# plt.tight_layout( )
# plt.grid( )
plt.show( )
# plt.savefig('../data/trajectory_for_different_launch_declination/trajectoryPlot_'+str(figureIndex+1)+'.png', dpi=300)

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""