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
# Connect to SQLite database.
try:
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/"
        #                             + "spherical_asteroid/"
        #                             + "longestEdge.db" )
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/"
        #                             + "multiple_launch_velocity_with_perturbations/"
        #                             + "simulation_time_9_months/"
        #                             + "3.2Density_1cmRadius/"
        #                             + "longestEdgePerturbations.db" )
        # database = sqlite3.connect("../data/VandV/"
        #                         + "new/"
        #                         + "launch_velocity/"
        #                         + "leadingEdge2.db")
        database = sqlite3.connect("../data/VandV/"
                                + "new/"
                                + "orbital_dynamics/"
                                + "leadingEdge3.db")
        # database = sqlite3.connect("../data/VandV/"
        #                         + "new/"
        #                         + "orbital_dynamics/"
        #                         + "sphericalAsteroid.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

data = pd.read_sql( "SELECT     position_x,                                                 \
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
                     WHERE      ROUND( launch_azimuth ) = 135.0                             \
                     AND        ROUND( launch_declination ) = 30.0                          \
                     AND        ROUND( initial_velocity_magnitude ) = 10.0;",               \
                     database )

data.columns = [ 'x',                                                   \
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

# upto 7 significant digits after the decimal point in jacobian and energy
for index in range( 0, len(jacobi) ):
    jacobi[index] = float( "{0:.7f}".format(jacobi[index]) )
    energy[index] = float( "{0:.7f}".format(energy[index]) )

## Set up the figure
fig = plt.figure( figsize=(6, 6) )
gs = gridspec.GridSpec( 2, 1 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
plt.suptitle("CDE asteroid\n"
             + "Launch velocity = " + str(vel_mag[ 0 ]) + " [m/s], "
             + "Azimuth = " + str(azimuth[ 0 ]) + " [deg], "
             + "Declination = " + str(declination[ 0 ]) + " [deg]")

## plot jacobian and the energy for a single trajectory
ax1.plot( t, jacobi )
ax1.grid(True)
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Jacobi Integral')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
# ax1.legend( ).draggable( )

# ax2.plot( t, trueEnergy )
ax2.plot( t, energy )
ax2.grid(True)
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Total Energy [$m^2/s^2$]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)

## Plot the ellipsoidal shape of the asteroid
# alpha = 20000.0
# beta = 7000.0
# gamma = 7000.0
# Wz = 0.00033118202125129593

# u = np.linspace(0, 2 * np.pi, 100)
# v = np.linspace(0, np.pi, 100)

# ellipsoid_x = alpha * np.outer(np.cos(u), np.sin(v))
# ellipsoid_y = beta * np.outer(np.sin(u), np.sin(v))
# ellipsoid_z = gamma * np.outer(np.ones(np.size(u)), np.cos(v))

# newColor = colors.cnames["slategray"]
# surf = ax1.plot_surface( ellipsoid_x, ellipsoid_y, ellipsoid_z,
#                          rstride=5, cstride=5 )
# # surf.set_facecolor( ( 0, 0, 1, 0.5 ) )
# surf.set_facecolor( newColor )
# surf.set_linewidth( 0.1 )

# ax1.hold( True )

# ## Plot 3D trajectory of the orbiting particle
# ax1.plot( x, y, z, zdir = 'z', color=colors.cnames["purple"] )

# ## indicate starting point
# ax1.text( x[0], y[0], z[0], 'start', size=10, zorder=1, color=colors.cnames["black"] )

# ## indicate ending point
# endIndex = np.size( x )
# ax1.text( x[endIndex-1], y[endIndex-1], z[endIndex-1],
#           'end', size=10, zorder=1,
#           color=colors.cnames["black"] )

# ## format axis and title
# ax1.set_xlabel('x [m]')
# ax1.set_ylabel('y [m]')
# ax1.set_zlabel('z [m]')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
# ax1.set_title( 'Particle trajectory around asteroid Eros (Body frame)' )

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
