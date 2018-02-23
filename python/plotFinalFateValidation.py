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

def plotEllipse( semiMajor, semiMinor, plotHandle ):
        angleRange = np.linspace(0, 2 * np.pi, 360)
        r = semiMajor * semiMinor / np.sqrt( ( semiMinor * np.cos( angleRange ) )**2
                                           + ( semiMajor * np.sin( angleRange ) )**2 )
        x = r * np.cos( angleRange )
        y = r * np.sin( angleRange )
        plotHandle.plot( x, y, label='Ellipsoid 2D projection' )
        plotHandle.grid( )

## constants
alpha = 20000.0
beta = 7000.0
gamma = 7000.0

# Start timer.
start_time = time.time( )

## Operations
# Connect to SQLite database.
try:
        # for escape and reimpact
        database = sqlite3.connect("../data/VandV/"
                                + "new/"
                                + "finalFateValidation.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

data = pd.read_sql( "SELECT     ROUND( trajectory_id ),                                     \
                                position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                velocity_x,                                                 \
                                velocity_y,                                                 \
                                velocity_z,                                                 \
                                ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth ),                                    \
                                ROUND( launch_declination ),                                \
                                time,                                                       \
                                sma,                                                        \
                                eccentricity,                                               \
                                inclination,                                                \
                                raan,                                                       \
                                aop,                                                        \
                                jacobi_integral,                                            \
                                total_energy                                                \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      escape_flag = 1;",                                          \
                     database )

data.columns = [ 'trajectory_id',                                       \
                 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'vx',                                                  \
                 'vy',                                                  \
                 'vz',                                                  \
                 'vel_mag',                                             \
                 'azimuth',                                             \
                 'declination',                                         \
                 'time',                                                \
                 'sma',                                                 \
                 'eccentricity',                                        \
                 'inclination',                                         \
                 'raan',                                                \
                 'aop',                                                 \
                 'jacobi_integral',                                     \
                 'total_energy' ]

trajectory_id   = data['trajectory_id']
x               = data[ 'x' ]
y               = data[ 'y' ]
z               = data[ 'z' ]
vx              = data[ 'vx' ]
vy              = data[ 'vy' ]
vz              = data[ 'vz' ]
vel_mag         = data['vel_mag']
azimuth         = data['azimuth']
declination     = data['declination']
t               = data[ 'time' ]
sma             = data['sma']
eccentricity    = data['eccentricity']
inclination     = data['inclination']
raan            = data['raan']
aop             = data['aop']
jacobi          = data[ 'jacobi_integral' ]
energy          = data['total_energy']

## Set up the figure
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
plt.suptitle("CDE asteroid, Launch Latitude = $60^o$, Launch Longitude = $30^o$")

for index in range( 0, len( trajectory_id ) ):
    eccentricityPlot = eccentricity[ index ]
    currentLaunchVelocity = vel_mag[ index ]
    currentLaunchAzimuth = azimuth[ index ]
    currentLaunchDeclination = declination[ index ]
    energyPlot = energy[ index ]

    ax1.plot( eccentricityPlot, energyPlot, marker='s',
              label=str(currentLaunchVelocity) + ", " + str(currentLaunchAzimuth) + ", " + str(currentLaunchDeclination) )

ax1.grid(True)
ax1.set_xlabel('Eccentricity')
ax1.set_ylabel('Keplerian Energy [$m^2/s^2$]')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax1.legend( title="Launch velocity, Azimuth, and Declination" ).draggable( )

## now plot the crash case validation
data = pd.read_sql( "SELECT     ROUND( trajectory_id ),                                     \
                                position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                velocity_x,                                                 \
                                velocity_y,                                                 \
                                velocity_z,                                                 \
                                ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth ),                                    \
                                ROUND( launch_declination ),                                \
                                time,                                                       \
                                sma,                                                        \
                                eccentricity,                                               \
                                inclination,                                                \
                                raan,                                                       \
                                aop,                                                        \
                                jacobi_integral,                                            \
                                total_energy                                                \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      crash_flag = 1                                              \
                     AND        ROUND( initial_velocity_magnitude ) = 13;",                 \
                     database )

data.columns = [ 'trajectory_id',                                       \
                 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'vx',                                                  \
                 'vy',                                                  \
                 'vz',                                                  \
                 'vel_mag',                                             \
                 'azimuth',                                             \
                 'declination',                                         \
                 'time',                                                \
                 'sma',                                                 \
                 'eccentricity',                                        \
                 'inclination',                                         \
                 'raan',                                                \
                 'aop',                                                 \
                 'jacobi_integral',                                     \
                 'total_energy' ]

trajectory_id   = data['trajectory_id']
x               = data[ 'x' ]
y               = data[ 'y' ]
z               = data[ 'z' ]
vx              = data[ 'vx' ]
vy              = data[ 'vy' ]
vz              = data[ 'vz' ]
vel_mag         = data['vel_mag']
azimuth         = data['azimuth']
declination     = data['declination']
t               = data[ 'time' ]
sma             = data['sma']
eccentricity    = data['eccentricity']
inclination     = data['inclination']
raan            = data['raan']
aop             = data['aop']
jacobi          = data[ 'jacobi_integral' ]
energy          = data['total_energy']

## Set up the figure
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
plt.suptitle("CDE asteroid, Launch Latitude = $60^o$, Launch Longitude = $30^o$")

for index in range( 0, len( trajectory_id ) ):
    currentLaunchVelocity = vel_mag[ index ]
    currentLaunchAzimuth = azimuth[ index ]
    currentLaunchDeclination = declination[ index ]

    cdeSolution = x[index]**2 / alpha**2 + y[index]**2 / beta**2 + z[index]**2 / gamma**2 - 1.0
    ax1.plot( t[index], cdeSolution, marker='s',
              label=str(currentLaunchVelocity) + ", " + str(currentLaunchAzimuth) + ", " + str(currentLaunchDeclination) )

    # cdeSolution = x**2 / alpha**2 + y**2 / beta**2 + z**2 / gamma**2 - 1.0
    # ax1.plot( t, cdeSolution, marker='s',
    #           label=str(currentLaunchVelocity) + ", " + str(currentLaunchAzimuth) + ", " + str(currentLaunchDeclination) )

ax1.grid(True)
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Solution to CDE equation')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax1.legend( title="Launch velocity, Azimuth, and Declination" ).draggable( )

## crash case plot 2
data = pd.read_sql( "SELECT     ROUND( trajectory_id )                                      \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      crash_flag = 1;",                                           \
                     database )

data.columns = ['trajectory_id']
trajectory_id = data['trajectory_id']
trajectory_id = tuple( trajectory_id.tolist( ) )

data = pd.read_sql( "SELECT     ROUND( trajectory_id ),                                     \
                                position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                velocity_x,                                                 \
                                velocity_y,                                                 \
                                velocity_z,                                                 \
                                ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth ),                                    \
                                ROUND( launch_declination ),                                \
                                time,                                                       \
                                sma,                                                        \
                                eccentricity,                                               \
                                inclination,                                                \
                                raan,                                                       \
                                aop,                                                        \
                                jacobi_integral,                                            \
                                total_energy                                                \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      trajectory_id IN " + str( trajectory_id ) + "               \
                     AND        ROUND( initial_velocity_magnitude ) = 13;",                 \
                     database )

data.columns = [ 'trajectory_id',                                       \
                 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'vx',                                                  \
                 'vy',                                                  \
                 'vz',                                                  \
                 'vel_mag',                                             \
                 'azimuth',                                             \
                 'declination',                                         \
                 'time',                                                \
                 'sma',                                                 \
                 'eccentricity',                                        \
                 'inclination',                                         \
                 'raan',                                                \
                 'aop',                                                 \
                 'jacobi_integral',                                     \
                 'total_energy' ]

trajectory_id   = data['trajectory_id']
x               = data[ 'x' ]
y               = data[ 'y' ]
z               = data[ 'z' ]
vx              = data[ 'vx' ]
vy              = data[ 'vy' ]
vz              = data[ 'vz' ]
vel_mag         = data['vel_mag']
azimuth         = data['azimuth']
declination     = data['declination']
t               = data[ 'time' ]
sma             = data['sma']
eccentricity    = data['eccentricity']
inclination     = data['inclination']
raan            = data['raan']
aop             = data['aop']
jacobi          = data[ 'jacobi_integral' ]
energy          = data['total_energy']

## Set up the figure
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
plt.suptitle("CDE asteroid, Launch Latitude = $60^o$, Launch Longitude = $30^o$")

uniqueTrajectoryIds = np.unique( trajectory_id )
for index in range( 0, len( uniqueTrajectoryIds ) ):
    currentTrajectoryId = uniqueTrajectoryIds[ index ]
    data_indices = np.where( trajectory_id == currentTrajectoryId )
    data_indices = data_indices[ 0 ]

    currentLaunchVelocity = (vel_mag[ data_indices ]).tolist( )
    currentLaunchVelocity = currentLaunchVelocity[ -1 ]

    currentLaunchAzimuth = azimuth[ data_indices ].tolist( )
    currentLaunchAzimuth = currentLaunchAzimuth[ -1 ]

    currentLaunchDeclination = declination[ data_indices ].tolist( )
    currentLaunchDeclination = currentLaunchDeclination[ -1 ]

    cdeSolution = x[data_indices]**2 / alpha**2 + y[data_indices]**2 / beta**2 + z[data_indices]**2 / gamma**2 - 1.0
    ax1.plot( t[data_indices], cdeSolution,
              label=str(currentLaunchVelocity) + ", " + str(currentLaunchAzimuth) + ", " + str(currentLaunchDeclination) )

ax1.grid(True)
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Solution to CDE equation')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax1.legend( title="Launch velocity, Azimuth, and Declination" ).draggable( )

if database:
    database.close( )

try:
        # for capture
        database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/"
                                    + "multiple_launch_velocity_with_perturbations/"
                                    + "simulation_time_9_months/"
                                    # + "3.2Density_1cmRadius/"
                                    + "test.db" )
except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

data = pd.read_sql( "SELECT     ROUND( trajectory_id ),                                     \
                                position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                inertial_position_x,                                        \
                                inertial_position_y,                                        \
                                inertial_position_z,                                        \
                                ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth ),                                    \
                                ROUND( launch_declination ),                                \
                                time,                                                       \
                                sma,                                                        \
                                eccentricity,                                               \
                                inclination,                                                \
                                raan,                                                       \
                                aop,                                                        \
                                jacobi_integral,                                            \
                                total_energy                                                \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ROUND( initial_velocity_magnitude ) = 10                    \
                     AND        ROUND( launch_azimuth ) = 45                                \
                     AND        ROUND( initial_solar_phase_angle ) = 315;",                 \
                     database )

data.columns = [ 'trajectory_id',                                       \
                 'x',                                                   \
                 'y',                                                   \
                 'z',                                                   \
                 'inertial_x',                                          \
                 'inertial_y',                                          \
                 'inertial_z',                                          \
                 'vel_mag',                                             \
                 'azimuth',                                             \
                 'declination',                                         \
                 'time',                                                \
                 'sma',                                                 \
                 'eccentricity',                                        \
                 'inclination',                                         \
                 'raan',                                                \
                 'aop',                                                 \
                 'jacobi_integral',                                     \
                 'total_energy' ]

trajectory_id   = data['trajectory_id']
x               = data[ 'x' ]
y               = data[ 'y' ]
z               = data[ 'z' ]
inertial_x      = data[ 'inertial_x' ]
inertial_y      = data[ 'inertial_y' ]
inertial_z      = data[ 'inertial_z' ]
vel_mag         = data['vel_mag']
azimuth         = data['azimuth']
declination     = data['declination']
t               = data[ 'time' ]
sma             = data['sma']
eccentricity    = data['eccentricity']
inclination     = data['inclination']
raan            = data['raan']
aop             = data['aop']
jacobi          = data[ 'jacobi_integral' ]
energy          = data['total_energy']

## set up the figure
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ], projection='3d' )

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

newColor = colors.cnames["brown"]
surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
                         rstride=5, cstride=5, alpha=0.8,
                         color=newColor )
surf.set_linewidth( 1.0 )

particleRange = np.sqrt( x**2 + y**2 + z**2 )
data_indices = np.where( particleRange <= 1.5 * alpha )
data_indices = data_indices[ 0 ]
xPlot = x[ data_indices ]
yPlot = y[ data_indices ]
zPlot = z[ data_indices ]
ax1.scatter( xPlot, yPlot, zPlot, alpha=0.5 )

ax1.grid(True)
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')
# ax1.set_xlim(-20000, 20000)
# ax1.set_ylim(-20000, 20000)
# ax1.set_zlim(-20000, 20000)
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)

## Set up the figure
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

data_indices = np.where( particleRange >= 10.0 * alpha )
data_indices = data_indices[ 0 ]

eccentricityPlot = eccentricity[ data_indices ]
energyPlot = energy[ data_indices ]

ax1.scatter( eccentricityPlot, energyPlot, s=1 )

ax1.grid(True)
ax1.set_xlabel('Eccentricity')
ax1.set_ylabel('Energy [$m^2/s^2$]')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)

## Set up the figure
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

ax1.plot(inertial_x, inertial_y)
ax1.grid(True)
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)

## Set up the figure
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

plotEllipse( alpha, beta, ax1 )
ax1.plot(x, y, lw=0.02)
ax1.grid(True)
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)

## Set up the figure
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
plt.suptitle("CDE asteroid, Launch Latitude = $0^o$, Launch Longitude = $0^o$")

rangeToRegolith = np.sqrt( x**2 + y**2 + z**2 )

newAlpha = alpha + 1.0e-2
newBeta = beta + 1.0e-2
newGamma = gamma + 1.0e-2
cdeSolution = x**2 / newAlpha**2 + y**2 / newBeta**2 + z**2 / newGamma**2 - 1.0

for index in range( 0, len( cdeSolution ) ):
    if cdeSolution[ index ] <= 0.0 and index > 1:
        print cdeSolution[index]
        print t[index]
        print "Error: Particle crashed in capture case plot"
        sys.exit( )

ax1.plot( t, eccentricity )
# ax1.plot( rangeToRegolith, eccentricity )
ax1.grid(True)
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Eccentricity')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)

fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
plt.suptitle("CDE asteroid, Launch Latitude = $0^o$, Launch Longitude = $0^o$")

ax1.plot( t, energy )
# ax1.plot( rangeToRegolith, energy )
ax1.grid(True)
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Energy [$m^2/s^2$]')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)

fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
plt.suptitle("CDE asteroid, Launch Latitude = $0^o$, Launch Longitude = $0^o$")

ax1.plot( t, cdeSolution )
# ax1.plot( rangeToRegolith, cdeSolution )
ax1.grid( True )
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Solution to CDE equation')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)

fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
plt.suptitle("CDE asteroid, Launch Latitude = $0^o$, Launch Longitude = $0^o$")

ax1.plot( t, rangeToRegolith )
# ax1.plot( rangeToRegolith, cdeSolution )
ax1.grid( True )
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Range [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)

# Stop timer
end_time = time.time( )

## Show the plot
plt.show( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
