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

## plots for the non perturbing case and the same initial conditions
# Connect to SQLite database.
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
                     WHERE      ROUND( launch_azimuth ) = 85.0                      \
                     AND        ROUND( initial_velocity_magnitude ) = 14.0;",       \
                     database )

data.columns = [ 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'velocity_magnitude',                                  \
                 'inertial_x',                                          \
                 'inertial_y',                                          \
                 'inertial_z',                                          \
                 'inertial_velocity_x',                                 \
                 'inertial_velocity_y',                                 \
                 'inertial_velocity_z',                                 \
                 'launch_azimuth',                                      \
                 'time' ]

x                   = data[ 'x' ]
y                   = data[ 'y' ]
z                   = data[ 'z' ]
velocityMagnitude   = data[ 'velocity_magnitude' ]
inertial_x          = data[ 'inertial_x' ]
inertial_y          = data[ 'inertial_y' ]
inertial_z          = data[ 'inertial_z' ]
inertial_xVel       = data[ 'inertial_velocity_x' ]
inertial_yVel       = data[ 'inertial_velocity_y' ]
inertial_zVel       = data[ 'inertial_velocity_z' ]
launchAzimuth       = data[ 'launch_azimuth' ]

t_noSolarPerturbations = data[ 'time' ]

if database:
    database.close( )

# get the end point for a data array
endIndex = np.size( x )

string1 = 'Particle trajectory projection around asteroid Eros (Body fixed frame) \n'
string2 = '$V_{initial}$=' + str( velocityMagnitude[ 0 ] ) + '[m/s], '
string3 = 'Launch azimuth=' + str( launchAzimuth[ 0 ] ) + '[deg], '
string4 = 'time=' + str( t_noSolarPerturbations[ endIndex-1 ] / (60.0*60.0) ) + '[hrs]'
topTitleBodyFrame = string1 + string2 + string3 + string4

fig = plt.figure( )
plt.suptitle( topTitleBodyFrame )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

trajectoryColor = colors.cnames["purple"]
textColor       = colors.cnames["black"]
startColor      = colors.cnames["darkgreen"]
endColor        = colors.cnames["darkred"]

ax1.plot( x, y, color=trajectoryColor )

## indicate starting point
ax1.text( x[0], y[0], 'start', size=12, color=startColor )

## indicate ending point
ax1.text( x[endIndex-1], y[endIndex-1], 'end', size=12, color=endColor )

## format axis and title
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
ax1.grid( True )
ax1.axis('equal')


## plot for the perturbing case
# Connect to SQLite database.
try:
    database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                               + "multiple_launch_velocity_with_perturbations/"
                               + "simulation_time_9_months/"
                               + "3.2Density_1cmSize/longestEdgePerturbations.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data for the perturbing case now...\n"

# manually set time for data extraction
lowerTimeDays = 0.0
lowerTime = lowerTimeDays * 24.0 * 60.0 * 60.0

# auto set time from the non-perturbing results
upperTime = t_noSolarPerturbations[ endIndex-1 ]
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
                                ROUND( time )                                                   \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( launch_azimuth ) = 85.0                                  \
                     AND        ROUND( initial_velocity_magnitude ) = 14.0                      \
                     AND        ROUND( initial_solar_phase_angle ) = 135.0;",                   \
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

## get data for the sun's trajectory
data2 = pd.read_csv("../data/regolith_launched_from_longest_edge/"
                    + "multiple_launch_velocity_with_perturbations/"
                    + "simulation_time_9_months/"
                    + "3.2Density_1cmSize/sunEphemeris_phase135.csv")

sun_xBody           = data2["x_body_frame"].values
sun_yBody           = data2["y_body_frame"].values
sun_zBody           = data2["z_body_frame"].values
sun_xInertial       = data2["x_inertial_frame"].values
sun_yInertial       = data2["y_inertial_frame"].values
sun_zInertial       = data2["z_inertial_frame"].values
sunLongitude        = data2["ta"].values
sunKeplerSolverTime = data2["t"].values
sunEndIndex         = len( sun_xBody ) - 1

print "Processing data now...\n"

print "Solar phase at start of the segment = " + str( changingSolarPhase[ 0 ] ) + " [deg]"

## sanity check for sun's trajectory from the CSV file (plot sun logitude)
fig = plt.figure( figsize=( 10, 10 ) )
gs = gridspec.GridSpec( 2, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

ax1.plot( ( sunKeplerSolverTime/( 24.0 * 60.0 * 60.0 ) ), sunLongitude )
ax1.grid(True)
ax1.set_ylabel('Solar longitude [deg]')
ax1.set_xlabel('Time [days]')
ax1.ticklabel_format( style='sci', axis='x', scilimits=(0,0), useOffset=False )

ax2.plot( sun_xBody, sun_yBody )
ax2.plot( sun_xBody[0], sun_yBody[0],
          marker='*', markersize=5, color='green', label='Sun start' )
ax2.plot( sun_xBody[sunEndIndex], sun_yBody[sunEndIndex],
          marker='*', markersize=5, color='black', label='Sun end' )
ax2.grid(True)
ax2.legend( ).draggable( )
ax2.set_ylabel('y [m]')
ax2.set_xlabel('x [m]')
ax2.set_title('Rotating frame')
ax2.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )

ax3.plot( sun_xInertial, sun_yInertial )
ax3.plot( sun_xInertial[0], sun_yInertial[0],
          marker='*', markersize=5, color='green', label='Sun start' )
ax3.plot( sun_xInertial[sunEndIndex], sun_yInertial[sunEndIndex],
          marker='*', markersize=5, color='black', label='Sun end' )
ax3.grid(True)
ax3.legend( ).draggable( )
ax3.set_ylabel('y [m]')
ax3.set_xlabel('x [m]')
ax3.set_title('Inertial frame')
ax3.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )

## Set up the figure for the actual plot
fig = plt.figure( figsize=( 10, 10 ) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
plt.suptitle( 'Ellipsoid longest edge, Body frame plot\n'
              + '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]\n'
              + 'Time = ' + str( lowerTimeDays ) + ' to ' + str( upperTimeDays ) + ' [days]',
              fontsize=10 )

### main segment of the code
# Draw the entire trajectory path for the regolith
regolithTrajectory = ax1.plot( position_x, position_y,
                               linewidth=0.1, color=colors.cnames['purple'],
                               label='Regolith trajectory' )

## data indices
xyRange = np.sqrt( position_x**2 + position_y**2 )
limitLowerTrajectory = np.sqrt( 0.5e5**2 + 0.5e5**2 )
limitUpperTrajectory = np.sqrt( 1.0e5**2 + 1.0e5**2 )

data_indices_1 = np.where( xyRange <= limitLowerTrajectory )
data_indices_1 = data_indices_1[ 0 ]
data_indices_1 = data_indices_1[ 0::500 ]

data_indices_2 = np.where( (xyRange > limitLowerTrajectory) & (xyRange < limitUpperTrajectory) )
data_indices_2 = data_indices_2[ 0 ]
data_indices_2 = data_indices_2[ 0::200 ]

data_indices_3 = np.where( xyRange > limitUpperTrajectory )
data_indices_3 = data_indices_3[ 0 ]

data_indices = data_indices_1.tolist( ) + data_indices_2.tolist( ) + data_indices_3.tolist( )

plotTimes = t[data_indices]
plotTimes = plotTimes.tolist( )

sunTimeValues = []
sunDataIndices = []
for index in range( 0, len( plotTimes ) ):
    # sunDataIndicesTemp = np.where( sunKeplerSolverTime == plotTimes[ index ] )
    # sunDataIndicesTemp = sunDataIndicesTemp[ 0 ]
    sunDataIndicesTemp = (np.abs(sunKeplerSolverTime - plotTimes[index])).argmin( )
    sunDataIndices.append( sunDataIndicesTemp )
    sunTimeValues.append( sunKeplerSolverTime[ sunDataIndicesTemp ] )

# sunTimeValues = np.concatenate( sunTimeValues, axis=0 )
# sunTimeValues = sunTimeValues.tolist( )

plotStepSize = 1

## srp vector plot
y = position_y[ data_indices ].tolist( )
x = position_x[ data_indices ].tolist( )
z = position_z[ data_indices ].tolist( )
srp_xPlot = srp_x[ data_indices ].tolist( )
srp_yPlot = srp_y[ data_indices ].tolist( )
srpVector = []
for index in range( 0, len( srp_xPlot ), plotStepSize ):
    currentVector = ax1.quiver( x[index], y[index],
                                srp_xPlot[index], srp_yPlot[index],
                                angles='xy',
                                width=0.002,
                                headwidth=4,
                                headlength=4,
                                headaxislength=4,
                                linewidths=0.5,
                                pivot='tail',
                                color='red',
                                label='SRP' if index == 0 else '_nolegend_' )

    srpVector.append( currentVector )

## solar third body effect
solarTide_xPlot = solarTide_x[ data_indices ].tolist( )
solarTide_yPlot = solarTide_y[ data_indices ].tolist( )
for index in range( 0, len( solarTide_xPlot ), plotStepSize ):
    currentVector = ax1.quiver( x[index], y[index],
                                solarTide_xPlot[index], solarTide_yPlot[index],
                                angles='xy',
                                width=0.002,
                                headwidth=4,
                                headlength=4,
                                headaxislength=4,
                                linewidths=0.5,
                                pivot='tail',
                                color='green',
                                label='Solar-TBE' if index == 0 else '_nolegend_' )

## sort out sun data to be plotted
# sunDataIndices = np.concatenate( sunDataIndices, axis=0 ).tolist( )
sun_xPlot = sun_xBody[ sunDataIndices ]
sun_yPlot = sun_yBody[ sunDataIndices ]
sun_zPlot = sun_zBody[ sunDataIndices ]

## scatter plot the sun and regolith's position corresponding to the same given time
colors = plt.cm.Vega20( np.linspace( 0, 1, len( x ) ) )
for index in range( 0, len( x ), plotStepSize ):
    ax1.scatter( x[index], y[index], s=40, c=colors[index],
                 marker='o', edgecolors='face',
                 label='Regolith' if index == 0 else '_nolegend_' )

ax1.grid(True)
ax1.set_xlabel( 'x [m]' )
ax1.set_ylabel( 'y [m]' )
# ax1.set_xlim( min( position_x ), max( position_x ) )
# ax1.set_ylim( min( position_y ), max( position_y ) )
ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
leg = ax1.legend( ).draggable( )
# ax1.set_aspect( 1 )

if database:
    database.close( )

## do the same analysis for inertial frame orbits
sun_xInertialPlot = sun_xInertial[ sunDataIndices ]
sun_yInertialPlot = sun_yInertial[ sunDataIndices ]
sun_zInertialPlot = sun_zInertial[ sunDataIndices ]

regolith_xInertialPlot = inertial_position_x[ data_indices ]
regolith_yInertialPlot = inertial_position_y[ data_indices ]
regolith_zInertialPlot = inertial_position_z[ data_indices ]

# compute SRP in the inertial frame
rho = 1.0
solarConstant = 1.0e17
muSun = 1.32712440018 * 1.0e+20

regolithGrainDensity = 3.2 * 1.0e3;
regolithGrainRadius = 1.0 * 1.0e-2;
regolithCrossSectionalArea = np.pi * regolithGrainRadius**2
regolithMass = (4.0 / 3.0) * np.pi * regolithGrainRadius**3 * regolithGrainDensity
areaToMassRatio = regolithCrossSectionalArea / regolithMass;

multiplicationConstant = -1.0 * ( 1.0 + rho ) * solarConstant * areaToMassRatio;

DminusR_xInertial = sun_xInertialPlot - regolith_xInertialPlot
DminusR_yInertial = sun_yInertialPlot - regolith_yInertialPlot
DminusR_zInertial = sun_zInertialPlot - regolith_zInertialPlot

DminusRcube = ( np.sqrt( DminusR_xInertial**2 + DminusR_yInertial**2 + DminusR_zInertial**2 ) )**3

srp_xInertial = multiplicationConstant * DminusR_xInertial / ( DminusRcube )
srp_yInertial = multiplicationConstant * DminusR_yInertial / ( DminusRcube )
srp_zInertial = multiplicationConstant * DminusR_zInertial / ( DminusRcube )

# compute STBE for the inertial frame
RminusD_xInertial = regolith_xInertialPlot - sun_xInertialPlot
RminusD_yInertial = regolith_yInertialPlot - sun_yInertialPlot
RminusD_zInertial = regolith_zInertialPlot - sun_zInertialPlot

RminusD_MagnitudeCubeInertial = ( np.sqrt( RminusD_xInertial**2 + RminusD_yInertial**2 + RminusD_zInertial**2 ) )**3

sun_MagnitudeCubeInertial = ( np.sqrt( sun_xInertialPlot**2 + sun_yInertialPlot**2 + sun_zInertialPlot**2 ) )**3

stbe_xInertial = -muSun * ( RminusD_xInertial / RminusD_MagnitudeCubeInertial + sun_xInertialPlot / sun_MagnitudeCubeInertial )
stbe_yInertial = -muSun * ( RminusD_yInertial / RminusD_MagnitudeCubeInertial + sun_yInertialPlot / sun_MagnitudeCubeInertial )
stbe_zInertial = -muSun * ( RminusD_zInertial / RminusD_MagnitudeCubeInertial + sun_zInertialPlot / sun_MagnitudeCubeInertial )

# set up the figure now
fig = plt.figure( figsize=(20, 20) )
gs = gridspec.GridSpec( 1, 2 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
plt.suptitle( 'Ellipsoid longest edge, Inertial frame plot\n'
              + '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]\n'
              + 'Time = ' + str( lowerTimeDays ) + ' to ' + str( upperTimeDays ) + ' [days]',
              fontsize=10 )

# plot the inertial frame trajectory for unperturbed case
ax1.plot( inertial_x, inertial_y,
          color=trajectoryColor, linewidth=0.5,
          label='Regolith trajectory' )

# indicate starting point
ax1.text( inertial_x[0], inertial_y[0], 'start', size=12, color=startColor )

# indicate ending point
ax1.text( inertial_x[endIndex-1], inertial_y[endIndex-1], 'end', size=12, color=endColor )

# format axis and title
ax1.set_title( 'No Solar perturbations' )
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
ax1.grid( True )
ax1.axis('equal')
ax1.legend( ).draggable( )

# plot inertial trajectory for perturbed case
ax2.plot( inertial_position_x, inertial_position_y,
          color=trajectoryColor, linewidth=0.5 )

# mark the regolith points and draw the srp vector
inertialSRPSanityCheckFlag = False
colors = plt.cm.Vega20( np.linspace( 0, 1, len( regolith_xInertialPlot ) ) )

regolith_xInertialPlot = regolith_xInertialPlot.tolist( )
regolith_yInertialPlot = regolith_yInertialPlot.tolist( )
srp_xInertial = srp_xInertial.tolist( )
srp_yInertial = srp_yInertial.tolist( )
sun_xInertialPlot = sun_xInertialPlot.tolist( )
sun_yInertialPlot = sun_yInertialPlot.tolist( )
stbe_xInertial = stbe_xInertial.tolist( )
stbe_yInertial = stbe_yInertial.tolist( )

for index in range( 0, len( regolith_xInertialPlot ), plotStepSize ):
    # ax2.scatter( regolith_xInertialPlot[index], regolith_yInertialPlot[index],
    #              s=40, c=colors[index],
    #              marker='o', edgecolors='face',
    #              label='Regolith' if index == 0 else '_nolegend_' )

    ax2.quiver( regolith_xInertialPlot[index], regolith_yInertialPlot[index],
                srp_xInertial[index], srp_yInertial[index],
                angles='xy',
                width=0.002,
                headwidth=4,
                headlength=4,
                headaxislength=4,
                linewidths=0.5,
                pivot='tail',
                color='red',
                label='SRP' if index == 0 else '_nolegend_' )

    # ax2.quiver( regolith_xInertialPlot[index], regolith_yInertialPlot[index],
    #             stbe_xInertial[index], stbe_yInertial[index],
    #             angles='xy',
    #             width=0.002,
    #             headwidth=4,
    #             headlength=4,
    #             headaxislength=4,
    #             linewidths=0.5,
    #             pivot='tail',
    #             color='green',
    #             label='STBE' if index == 0 else '_nolegend_' )

    if inertialSRPSanityCheckFlag == True:
        ax2.quiver( regolith_xInertialPlot[index], regolith_yInertialPlot[index],
                    sun_xInertialPlot[index], sun_yInertialPlot[index],
                    angles='xy',
                    width=0.002,
                    headwidth=4,
                    headlength=4,
                    headaxislength=4,
                    linewidths=0.5,
                    pivot='tail',
                    color='black',
                    label="Sun's direction" if index == 0 else "_nolegend_" )

ax2.grid(True)
ax2.set_title( 'With Solar perturbations' )
ax2.set_xlabel( 'x [m]' )
ax2.set_ylabel( 'y [m]' )
ax2.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
ax2.axis('equal')
ax2.legend( ).draggable( )

# plot inertial velocity comparison
fig = plt.figure( )
ax1 = plt.subplot( 111 )
plt.suptitle( 'Ellipsoid longest edge, Inertial frame plot\n'
              + '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]\n'
              + 'Time = ' + str( lowerTimeDays ) + ' to ' + str( upperTimeDays ) + ' [days]',
              fontsize=10 )

ax1.plot( t_noSolarPerturbations / ( 24.0 * 60.0 * 60.0 ),
          np.sqrt( inertial_xVel**2 + inertial_yVel**2 + inertial_zVel**2 ),
          color='orange',
          linestyle='dotted',
          label='No solar perturbations' )

ax1.plot( t / ( 24.0 * 60.0 * 60.0 ),
          np.sqrt( inertial_velocity_x**2 + inertial_velocity_y**2 + inertial_velocity_z**2 ),
          color='cyan',
          label='SRP + STBE' )

ax1.grid(True)
ax1.set_xlabel( "Time [days]" )
ax1.set_ylabel( "Inertial velocity [m/s]" )
ax1.legend( ).draggable( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## show the plot
# kwargDict = dict( kw_args = { 'bbox_inches' : 'tight' } )
# anim.save( savePath, savefig_kwargs=kwargDict, bitrate=100 )
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
