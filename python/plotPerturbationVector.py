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

def plotEllipse( semiMajor, semiMinor, angleRange, plotHandle ):
        r = semiMajor * semiMinor / np.sqrt( ( semiMinor * np.cos( angleRange ) )**2
                                           + ( semiMajor * np.sin( angleRange ) )**2 )
        x = r * np.cos( angleRange )
        y = r * np.sin( angleRange )
        plotHandle.plot( x, y )
        plotHandle.grid( )

## Operations
alpha = 10000.0
beta = 10000.0
gamma = 10000.0
Wz = 0.00033118202125129593
muAsteroid = 876514
muSun = 1.32712440018 * 1.0e+20
AU = 149597870700.0

## Operations
# Connect to SQLite database.
try:
    # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db" )
    # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )
    database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid_with_perturbations/longestEdge.db" )

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

## get data for the trajectory and the perturbations
data1 = pd.read_sql( "SELECT    solar_phase_angle,                                          \
                                position_x,                                                 \
                                position_y,                                                 \
                                position_z,                                                 \
                                srp_x,                                                      \
                                srp_y,                                                      \
                                srp_z,                                                      \
                                solarTide_x,                                                \
                                solarTide_y,                                                \
                                solarTide_z                                                 \
                     FROM       regolith_trajectory_results;",                              \
                     database )

data1.columns = [ 'solar_phase_angle',                                                      \
                  'pos_x',                                                                  \
                  'pos_y',                                                                  \
                  'pos_z',                                                                  \
                  'srp_x',                                                                  \
                  'srp_y',                                                                  \
                  'srp_z',                                                                  \
                  'solarTide_x',                                                            \
                  'solarTide_y',                                                            \
                  'solarTide_z' ]

solar_phase_angle       = data1[ 'solar_phase_angle' ]
pos_x                   = data1[ 'pos_x' ]
pos_y                   = data1[ 'pos_y' ]
pos_z                   = data1[ 'pos_z' ]
srp_x                   = data1[ 'srp_x' ]
srp_y                   = data1[ 'srp_y' ]
srp_z                   = data1[ 'srp_z' ]
solarTide_x             = data1[ 'solarTide_x' ]
solarTide_y             = data1[ 'solarTide_y' ]
solarTide_z             = data1[ 'solarTide_z' ]

## get data for the sun's trajectory
data = pd.read_csv( '../data/regolith_launched_from_longest_edge/spherical_asteroid_with_perturbations/launch_10ms_solarPhase_0Deg/sunEphemeris.csv' )

sun_xBody           = data["x_body_frame"].values
sun_yBody           = data["y_body_frame"].values
sun_zBody           = data["z_body_frame"].values
sunEndIndex = len( sun_xBody ) - 1

## plot trajectory and perturbing acceleration vectors
fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )
# ax2 = plt.subplot( gs[ 1 ] )
## set global plot title
plt.suptitle( 'Perturbation direction \n Spherical asteroid case, Body Frame, '
              + 'Solar Phase angle = ' + str( solar_phase_angle[ 0 ] ) + ' [deg]' )

theta = np.linspace(0, 2 * np.pi, 300)
plotEllipse( alpha, beta, theta, ax1 )

ax1.plot( pos_x, pos_y, c=colors.cnames['purple'], label='Regolith trajectory' )

indices = range( 0, len( pos_x ), 100 )
plot_pos_x = pos_x[ indices ]
plot_pos_y = pos_y[ indices ]
plot_pos_z = pos_z[ indices ]
plot_srp_x = srp_x[ indices ]
plot_srp_y = srp_y[ indices ]
plot_srp_z = srp_z[ indices ]
plot_solarTide_x = solarTide_x[ indices ]
plot_solarTide_y = solarTide_y[ indices ]
plot_solarTide_z = solarTide_z[ indices ]

ax1.quiver( plot_pos_x, plot_pos_y,
            plot_srp_x, plot_srp_y,
            pivot='tail',
            color=colors.cnames[ 'red' ], linestyles='solid', lw=1, label='SRP' )

ax1.quiver( plot_pos_x, plot_pos_y,
            plot_solarTide_x, plot_solarTide_y,
            pivot='tail',
            color=colors.cnames[ 'orange' ], linestyles='solid', lw=1, label='Sun Third-Body effect' )

ax1.quiver( 0.0, 0.0,
            sun_xBody[ 0 ], sun_yBody[ 0 ],
            pivot='tail',
            color=colors.cnames[ 'black' ], linestyles='solid', lw=1, label='Sun start position' )

ax1.quiver( 0.0, 0.0,
            sun_xBody[ sunEndIndex ], sun_yBody[ sunEndIndex ],
            pivot='tail',
            color=colors.cnames[ 'slategray' ], linestyles='solid', lw=1, label='Sun end position' )

ax1.plot( sun_xBody/(4.7e6), sun_yBody/(4.7e6), label='Scaled down Sun trajectory' )

ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.legend( ).draggable( )
ax1.grid(True)
ax1.axis('equal')

## plot sun's trajectory around the asteroid
# ax2.plot( sun_xBody, sun_yBody, color=colors.cnames['purple'] )

# ax2.scatter( sun_xBody[ 0 ], sun_yBody[ 0 ], s=200,                                            \
#              c=colors.cnames['orange'], marker=(5,1), label='Sun' )

# ax2.scatter( 0.0, 0.0, s=100, c=colors.cnames['black'], marker=(5,0), label='Eros' )
# ax2.set_xlabel('x [m]')
# ax2.set_ylabel('y [m]')
# ax2.set_title('Sun trajectory around the asteroid (Body frame)')
# ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
# ax2.legend( ).draggable( )
# ax2.grid(True)
# ax2.axis('equal')

##print out the values for perturbation in a table
fig = plt.figure( )
plt.suptitle('Variation in magnitude of Perturbing Accelerations')
ax3 = plt.subplot(211)
ax4 = plt.subplot(212)

srpMag          = np.sqrt( srp_x**2 + srp_y**2 + srp_z**2 )
solarTideMag    = np.sqrt( solarTide_x**2 + solarTide_y**2 + solarTide_z**2 )
radialDistance  = np.sqrt( pos_x**2 + pos_y**2 + pos_z**2 )

ax3.plot( radialDistance, srpMag, label='SRP magnitude' )
ax4.plot( radialDistance, solarTideMag, label='Sun Third-Body effect magnitude' )

ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax3.set_xlabel('Radial distance to regolith from asteroid centre [m]')
ax3.set_ylabel( 'Perturbing acceleration magnitude $[m^2/s^2]$' )
# ax3.set_yscale('log')
ax3.legend( ).draggable( )
ax3.grid(True)

ax4.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax4.set_xlabel('Radial distance to regolith from asteroid centre [m]')
ax4.set_ylabel( 'Perturbing acceleration magnitude $[m^2/s^2]$' )
# ax4.set_yscale('log')
ax4.legend( ).draggable( )
ax4.grid(True)

## close the database
if database:
    database.close( )

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
