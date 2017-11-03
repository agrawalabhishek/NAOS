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
from operator import add

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
muSun = 1.32712440018 * 1.0e+20

trajectoryColor = colors.cnames["purple"]
textColor       = colors.cnames["black"]
startColor      = colors.cnames["darkgreen"]
endColor        = colors.cnames["darkred"]

## plot for the perturbing case
# Connect to SQLite database.
try:
    # database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
    #                            + "multiple_launch_velocity_with_perturbations/"
    #                            + "simulation_time_9_months/"
    #                            + "3.2Density_1cmSize/longestEdgePerturbations.db")
    database = sqlite3.connect("../data/regolith_launched_from_longest_edge/"
                               + "multiple_launch_velocity_with_perturbations/"
                               + "simulation_time_9_months/"
                               + "3.2Density_1cmSize/longestEdgePerturbations_singleCase_noSTBE.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

print "Extracting data for the perturbing case now...\n"

# manually set time for data extraction
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

## get data for the sun's trajectory
data2 = pd.read_csv("../data/regolith_launched_from_longest_edge/"
                    + "multiple_launch_velocity_with_perturbations/"
                    + "simulation_time_9_months/"
                    + "3.2Density_1cmSize/sunEphemeris_phase315.csv")

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

## sanity check for sun's trajectory from the CSV file (plot sun longitude)
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

### main segment of the code

## data indices
print "\nSelecting data indices for regolith's simulation data...\n"

## data indices based on distance of regolith from the asteroid
# xyRange = np.sqrt( inertial_position_x**2 + inertial_position_y**2 )
# limitLowerTrajectory = np.sqrt( 0.0e6**2 + 0.3e6**2 )
# limitUpperTrajectory = np.sqrt( 0.5e6**2 + 0.0e6**2 )

# data_indices_1 = np.where( xyRange <= limitLowerTrajectory )
# data_indices_1 = data_indices_1[ 0 ]
# data_indices_1 = data_indices_1[ 0::500 ]

# data_indices_2 = np.where( (xyRange > limitLowerTrajectory) & (xyRange < limitUpperTrajectory) )
# data_indices_2 = data_indices_2[ 0 ]
# data_indices_2 = data_indices_2[ 0::200 ]

# data_indices_3 = np.where( xyRange > limitUpperTrajectory )
# data_indices_3 = data_indices_3[ 0 ]

# data_indices = data_indices_1.tolist( ) + data_indices_2.tolist( ) + data_indices_3.tolist( )

## data indices based on magnitude of grav. acceleration acting on particle
gravMagnitudeLimit = 1.0e-6
gravMag = np.sqrt( gravAcc_x**2 + gravAcc_y**2 + gravAcc_z**2 )
data_indices = np.where( gravMag < gravMagnitudeLimit )
data_indices = data_indices[ 0 ]
data_indices = data_indices[ 0::100 ]

plotTimes = t[data_indices]
plotTimes = plotTimes.tolist( )
print "Data indices for regolith's motion selected...\n"

print "Selecting data indices for sun's motion around the asteroid...\n"
sunTimeValues = []
sunDataIndices = []
for index in range( 0, len( plotTimes ) ):
    # sunDataIndicesTemp = np.where( sunKeplerSolverTime == plotTimes[ index ] )
    # sunDataIndicesTemp = sunDataIndicesTemp[ 0 ]
    # sunDataIndices.append( sunDataIndicesTemp )
    # sunTimeValues.append( sunKeplerSolverTime[ sunDataIndicesTemp ] )
    sunDataIndicesTemp = (np.abs(sunKeplerSolverTime - plotTimes[index])).argmin( )
    sunDataIndices.append( sunDataIndicesTemp )
    sunTimeValues.append( sunKeplerSolverTime[ sunDataIndicesTemp ] )

print "All data indices obtained...\n"

# sunTimeValues = np.concatenate( sunTimeValues, axis=0 )
# sunTimeValues = sunTimeValues.tolist( )

timeDifference = [ np.abs( i - j ) for i, j in zip( plotTimes, sunTimeValues ) ]
for index in range( 0, len( timeDifference ) ):
    # over a minute of a difference
    if timeDifference[ index ] > 60.0:
        print timeDifference[ index ]
        print "\n Error: Large difference between sun and regolith time data \n"
        sys.exit( )

# step size for data points to be plotted

print "Plotting data now...\n"

if database:
    database.close( )

## plot 2d trajectory projections, in case it is needed for analysis
fig = plt.figure( figsize=( 7, 10 ) )
gs = gridspec.GridSpec( 3, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )

upperTimeDays = t[ len(t) - 1 ] / ( 24.0 * 60.0 * 60.0 )
plt.suptitle( 'Ellipsoid longest edge, Inertial frame trajectory\n'
              + '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]\n'
              + 'Time = ' + str( lowerTimeDays ) + ' to ' + str( upperTimeDays ) + ' [days]',
              fontsize=10 )

ax1.plot( inertial_position_x, inertial_position_y, color=trajectoryColor, linewidth=1.0 )
ax1.text( inertial_position_x[ 0 ],
          inertial_position_y[ 0 ],
          'start',
          color=startColor )
ax1.text( inertial_position_x[ len( inertial_position_x ) - 1 ],
          inertial_position_y[ len( inertial_position_y ) - 1 ],
          'end',
          color=endColor )
ax1.grid(True)
ax1.set_xlabel( 'x [m]' )
ax1.set_ylabel( 'y [m]' )
ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )

ax2.plot( inertial_position_y, inertial_position_z, color=trajectoryColor, linewidth=1.0 )
ax2.text( inertial_position_y[ 0 ],
          inertial_position_z[ 0 ],
          'start',
          color=startColor )
ax2.text( inertial_position_y[ len( inertial_position_x ) - 1 ],
          inertial_position_z[ len( inertial_position_y ) - 1 ],
          'end',
          color=endColor )
ax2.grid(True)
ax2.set_xlabel( 'y [m]' )
ax2.set_ylabel( 'z [m]' )
ax2.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )

ax3.plot( inertial_position_x, inertial_position_z, color=trajectoryColor, linewidth=1.0 )
ax3.text( inertial_position_x[ 0 ],
          inertial_position_z[ 0 ],
          'start',
          color=startColor )
ax3.text( inertial_position_x[ len( inertial_position_x ) - 1 ],
          inertial_position_z[ len( inertial_position_y ) - 1 ],
          'end',
          color=endColor )
ax3.grid(True)
ax3.set_xlabel( 'x [m]' )
ax3.set_ylabel( 'z [m]' )
ax3.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )

## analysis for inertial frame orbits
sun_xInertialPlot = sun_xInertial[ sunDataIndices ]
sun_yInertialPlot = sun_yInertial[ sunDataIndices ]
sun_zInertialPlot = sun_zInertial[ sunDataIndices ]

regolith_xInertialPlot = inertial_position_x[ data_indices ]
regolith_yInertialPlot = inertial_position_y[ data_indices ]
regolith_zInertialPlot = inertial_position_z[ data_indices ]

# compute SRP in the inertial frame
rho = 1.0
solarConstant = 1.0e17

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

## plot just the srp vector
regolith_xInertialPlot = regolith_xInertialPlot.tolist( )
regolith_yInertialPlot = regolith_yInertialPlot.tolist( )
regolith_zInertialPlot = regolith_zInertialPlot.tolist( )

srp_xInertial = srp_xInertial.tolist( )
srp_yInertial = srp_yInertial.tolist( )
srp_zInertial = srp_zInertial.tolist( )

sun_xInertialPlot = sun_xInertialPlot.tolist( )
sun_yInertialPlot = sun_yInertialPlot.tolist( )

fig = plt.figure( figsize=( 6, 6 ) )
ax1 = plt.subplot( 111 )

ax1.plot( inertial_position_x, inertial_position_y,
          color=trajectoryColor, linewidth=0.5, label='Regolith trajectory' )

ax1.text( inertial_position_x[ 0 ], inertial_position_y[ 0 ], 'start', color='green', fontsize=12 )
ax1.text( inertial_position_x[ len( inertial_position_x ) - 1 ],
          inertial_position_y[ len( inertial_position_x ) - 1 ],
          'end', color='red', fontsize=12 )

ax1.quiver( regolith_xInertialPlot, regolith_yInertialPlot,
            srp_xInertial, srp_yInertial,
            angles='xy',
            width=0.002,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color='red',
            label='SRP' )

ax1.grid(True)
# ax1.set_title( 'With Solar perturbations' )
ax1.set_xlabel( 'x [m]' )
ax1.set_ylabel( 'y [m]' )
ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
# ax1.axis('equal')
ax1.legend( ).draggable( )

upperTimeDays = t[ len(t) - 1 ] / ( 24.0 * 60.0 * 60.0 )
plt.suptitle( 'Ellipsoid longest edge, Inertial frame plot\n'
              + '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]\n'
              + 'Time = ' + str( lowerTimeDays ) + ' to ' + str( upperTimeDays ) + ' [days]',
              fontsize=10 )

## plot point mass gravity acceleration in inertial XY plot along with total perturbation vector
fig = plt.figure( figsize=( 6, 6 ) )
ax1 = plt.subplot( 111 )

# point mass gravity calculations
farAwayRangeMagnitude = np.sqrt( np.square( regolith_xInertialPlot ) \
                               + np.square( regolith_yInertialPlot ) \
                               + np.square( regolith_zInertialPlot ) )
farAwayRangeMagnitudeCube = farAwayRangeMagnitude**3

pointMass_xGravAcc = []
pointMass_yGravAcc = []
pointMass_zGravAcc = []
for index in range( 0, len( farAwayRangeMagnitudeCube ) ):
    pointMass_xGravAcc.append( ( -1.0 * mu * regolith_xInertialPlot[ index ] ) \
                                 / ( farAwayRangeMagnitudeCube[ index ] ) )

    pointMass_yGravAcc.append( ( -1.0 * mu * regolith_yInertialPlot[ index ] ) \
                                 / ( farAwayRangeMagnitudeCube[ index ] ) )

    pointMass_zGravAcc.append( ( -1.0 * mu * regolith_zInertialPlot[ index ] ) \
                                 / ( farAwayRangeMagnitudeCube[ index ] ) )

ax1.plot( inertial_position_x, inertial_position_y,
          color=trajectoryColor, linewidth=0.5, label='Regolith trajectory' )

ax1.text( inertial_position_x[ 0 ], inertial_position_y[ 0 ], 'start', color='green', fontsize=12 )
ax1.text( inertial_position_x[ len( inertial_position_x ) - 1 ],
          inertial_position_y[ len( inertial_position_x ) - 1 ],
          'end', color='red', fontsize=12 )

ax1.quiver( regolith_xInertialPlot, regolith_yInertialPlot,
            srp_xInertial, srp_yInertial,
            angles='xy',
            width=0.002,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color='red',
            label='SRP' )

ax1.quiver( regolith_xInertialPlot, regolith_yInertialPlot,
            pointMass_xGravAcc, pointMass_yGravAcc,
            angles='xy',
            width=0.002,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color='orange',
            label='gravity' )

ax1.grid( True )
ax1.set_xlabel( 'x [m]' )
ax1.set_ylabel( 'y [m]' )
ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
ax1.legend( ).draggable( )
upperTimeDays = t[ len(t) - 1 ] / ( 24.0 * 60.0 * 60.0 )
plt.suptitle( 'Ellipsoid longest edge, Inertial frame plot\n'
              + '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]\n'
              + 'Time = ' + str( lowerTimeDays ) + ' to ' + str( upperTimeDays ) + ' [days]',
              fontsize=10 )

## plot just the gravity vector for when the particle is far away from the asteroid
fig = plt.figure( figsize=( 6, 6 ) )
ax1 = plt.subplot( 111 )

ax1.plot( inertial_position_x, inertial_position_y,
          color=trajectoryColor, linewidth=0.5, label='Regolith trajectory' )

ax1.text( inertial_position_x[ 0 ], inertial_position_y[ 0 ], 'start', color='green', fontsize=12 )
ax1.text( inertial_position_x[ len( inertial_position_x ) - 1 ],
          inertial_position_y[ len( inertial_position_x ) - 1 ],
          'end', color='red', fontsize=12 )

ax1.quiver( regolith_xInertialPlot, regolith_yInertialPlot,
            pointMass_xGravAcc, pointMass_yGravAcc,
            angles='xy',
            width=0.002,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color='orange',
            label='gravity' )

ax1.grid( True )
ax1.set_xlabel( 'x [m]' )
ax1.set_ylabel( 'y [m]' )
ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
ax1.legend( ).draggable( )
upperTimeDays = t[ len(t) - 1 ] / ( 24.0 * 60.0 * 60.0 )
plt.suptitle( 'Ellipsoid longest edge, Inertial frame plot\n'
              + '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]\n'
              + 'Time = ' + str( lowerTimeDays ) + ' to ' + str( upperTimeDays ) + ' [days]',
              fontsize=10 )

## plot the net acceleration i.e. gravity plus all the perturbations
fig = plt.figure( figsize=( 6, 6 ) )
ax1 = plt.subplot( 111 )

ax1.plot( inertial_position_x, inertial_position_y,
          color=trajectoryColor, linewidth=0.5, label='Regolith trajectory' )

ax1.text( inertial_position_x[ 0 ], inertial_position_y[ 0 ], 'start', color='green', fontsize=12 )
ax1.text( inertial_position_x[ len( inertial_position_x ) - 1 ],
          inertial_position_y[ len( inertial_position_x ) - 1 ],
          'end', color='red', fontsize=12 )

totalAcceleration_xInertial = map( add, pointMass_xGravAcc, srp_xInertial )
totalAcceleration_yInertial = map( add, pointMass_yGravAcc, srp_yInertial )
totalAcceleration_zInertial = map( add, pointMass_zGravAcc, srp_zInertial )

ax1.quiver( regolith_xInertialPlot, regolith_yInertialPlot,
            totalAcceleration_xInertial, totalAcceleration_yInertial,
            angles='xy',
            width=0.002,
            headwidth=4,
            headlength=4,
            headaxislength=4,
            linewidths=0.5,
            pivot='tail',
            color=colors.cnames['royalblue'],
            label='Total acceleration' )

ax1.grid( True )
ax1.set_xlabel( 'x [m]' )
ax1.set_ylabel( 'y [m]' )
ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
ax1.legend( ).draggable( )
upperTimeDays = t[ len(t) - 1 ] / ( 24.0 * 60.0 * 60.0 )
plt.suptitle( 'Ellipsoid longest edge, Inertial frame plot\n'
              + '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]\n'
              + 'Time = ' + str( lowerTimeDays ) + ' to ' + str( upperTimeDays ) + ' [days]',
              fontsize=10 )

# sanity check
fig = plt.figure( figsize=( 6, 6 ) )
ax1 = plt.subplot( 111 )

timeToPlot = t[ data_indices ] / ( 24.0 * 60.0 * 60.0 )

totalPerturbingMagnitude = np.sqrt( np.square( srp_xInertial ) + \
                                    np.square( srp_yInertial ) + \
                                    np.square( srp_zInertial ) )

gravityMagnitude = np.sqrt( np.square( pointMass_xGravAcc ) + \
                            np.square( pointMass_yGravAcc ) + \
                            np.square( pointMass_zGravAcc ) )

totalAccelerationMagnitude = np.sqrt( np.square( totalAcceleration_xInertial ) + \
                                      np.square( totalAcceleration_yInertial ) + \
                                      np.square( totalAcceleration_zInertial ) )

ax1.scatter( timeToPlot, totalPerturbingMagnitude, color='red', label='SRP' )
ax1.scatter( timeToPlot, gravityMagnitude, color='orange', label='Gravitational acceleration' )
ax1.scatter( timeToPlot, totalAccelerationMagnitude, color=colors.cnames['royalblue'], label='Net acceleration' )

ax1.grid( True )
ax1.set_xlabel( 'Time [days]' )
ax1.set_ylabel( 'Magnitude [$m / s^2$]' )
ax1.set_yscale( 'log' )
# ax1.ticklabel_format( style='sci', axis='both', scilimits=(0,0), useOffset=False )
ax1.legend( ).draggable( )
upperTimeDays = t[ len(t) - 1 ] / ( 24.0 * 60.0 * 60.0 )
plt.suptitle( 'Ellipsoid longest edge, Inertial frame plot\n'
              + '$V_{launch}$ = ' + str( initialVelocityMagnitude[ 0 ] ) + ' [m/s], '
              + 'Launch azimuth = ' + str( launchAzimuth[ 0 ] ) + ' [deg], '
              + 'Solar phase = ' + str( solarPhase[ 0 ] ) + ' [deg]\n'
              + 'Time = ' + str( lowerTimeDays ) + ' to ' + str( upperTimeDays ) + ' [days]',
              fontsize=10 )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

## show the plot
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
