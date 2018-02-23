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

def plotCircle( radius, plotHandle ):
        angleRange = np.linspace( 0.0, 2.0*np.pi, 360.0 )
        x = radius * np.cos( angleRange )
        y = radius * np.sin( angleRange )
        plotHandle.plot( x, y )

## Operations
# Connect to SQLite database.
try:
    database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/"
                                + "multiple_launch_velocity/"
                                + "simulation_time_9_months/longestEdge.db" )

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

fig = plt.figure( figsize=(20, 20) )
gs = gridspec.GridSpec( 2, 3 )
ax1 = plt.subplot( gs[ 0 ] )
ax2 = plt.subplot( gs[ 1 ] )
ax3 = plt.subplot( gs[ 2 ] )
ax4 = plt.subplot( gs[ 3 ] )
ax5 = plt.subplot( gs[ 4 ] )
# ax6 = plt.subplot( gs[ 5 ] )
# ax7 = plt.subplot( gs[ 6 ] )

data2 = pd.read_sql( "SELECT    eccentricity,                                                   \
                                sma,                                                            \
                                inclination,                                                    \
                                raan,                                                           \
                                aop,                                                            \
                                ROUND( crash_flag ),                                            \
                                ROUND( escape_flag ),                                           \
                                ROUND( start_flag ),                                            \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth )                                         \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( initial_velocity_magnitude ) = 8.0;",                    \
                     database )

data2.columns = [ 'eccentricity',                                                       \
                  'sma',                                                                \
                  'inclination',                                                        \
                  'raan',                                                               \
                  'aop',                                                                \
                  'crash_flag',                                                         \
                  'escape_flag',                                                        \
                  'start_flag',                                                         \
                  'vel_mag',                                                            \
                  'launch_azimuth' ]

eccentricity                      =  data2[ 'eccentricity' ]
sma                               =  data2[ 'sma' ]
inclination                       =  data2[ 'inclination' ]
raan                              =  data2[ 'raan' ]
aop                               =  data2[ 'aop' ]
crashFlag                         =  data2[ 'crash_flag' ]
escapeFlag                        =  data2[ 'escape_flag' ]
startFlag                         =  data2[ 'start_flag' ]
current_velocity                  =  data2[ 'vel_mag' ]
launch_azimuth                    =  data2[ 'launch_azimuth' ]
current_velocity = current_velocity.tolist( )

# Close the database
if database:
    database.close( )

firstCrash = True
firstEscape = True

## Plot Ecc and SMA
uniqueAzimuths = np.unique( launch_azimuth )
for index in range( 0, len( uniqueAzimuths ) ):
    currentAzimuth = uniqueAzimuths[ index ]
    data_indices = np.where( launch_azimuth == currentAzimuth )
    data_indices = data_indices[ 0 ]

    current_crashFlag = crashFlag[ data_indices ]
    current_crashFlag = current_crashFlag.tolist( )
    current_escapeFlag = escapeFlag[ data_indices ]
    current_escapeFlag  = current_escapeFlag.tolist( )

    currentSMA = sma[ data_indices ]
    currentEcc = eccentricity[ data_indices ]
    currentInc = inclination[ data_indices ]
    currentRAAN = raan[ data_indices ]
    currentAOP = aop[ data_indices ]

    e_cos_w = currentEcc * np.cos( currentAOP * np.pi / 180.0 )
    e_cos_w = e_cos_w.tolist( )
    e_sin_w = currentEcc * np.sin( currentAOP * np.pi / 180.0 )
    e_sin_w = e_sin_w.tolist( )

    i_cos_raan = currentInc * np.cos( currentRAAN * np.pi / 180.0 )
    i_cos_raan = i_cos_raan.tolist( )
    i_sin_raan = currentInc * np.sin( currentRAAN * np.pi / 180.0 )
    i_sin_raan = i_sin_raan.tolist( )

    currentSMA = currentSMA.tolist( )
    currentEcc = currentEcc.tolist( )
    currentInc = currentInc.tolist( )
    currentRAAN = currentRAAN.tolist( )
    currentAOP = currentAOP.tolist( )

    if current_crashFlag[ -1 ] == 1.0:
        # particle with current azimuth had re-impacted - red
        ax1.scatter( currentAzimuth, currentSMA[0], s=5, c='red',                                   \
                     edgecolors='face',
                     label="Re-impact" if firstCrash else "__nolegend__" )
        firstCrash = False

        ax2.scatter( currentAzimuth, currentEcc[0], s=5, c='red',                                   \
                     edgecolors='face' )

        ax3.scatter( currentAzimuth, currentInc[0], s=5, c='red',                                   \
                     edgecolors='face' )

        ax4.scatter( currentAzimuth, currentRAAN[0], s=5, c='red',                                  \
                     edgecolors='face' )

        ax5.scatter( currentAzimuth, currentAOP[0], s=5, c='red',                                   \
                     edgecolors='face' )

        # ax6.scatter( e_cos_w[0], e_sin_w[0], s=5, c='red',                                          \
        #              edgecolors='face' )

        # ax7.scatter( i_cos_raan[0], i_sin_raan[0], s=5, c='red',                                    \
        #              edgecolors='face' )

    elif current_escapeFlag[ -1 ] == 1.0:
        # particle with current azimuth escaped - blue
        ax1.scatter( currentAzimuth, currentSMA[0], s=5, c='blue',                                  \
                     edgecolors='face',
                     label="Escape" if firstEscape else "__nolegend__" )
        firstEscape = False

        ax2.scatter( currentAzimuth, currentEcc[0], s=5, c='blue',                                  \
                     edgecolors='face' )

        ax3.scatter( currentAzimuth, currentInc[0], s=5, c='blue',                                  \
                     edgecolors='face' )

        ax4.scatter( currentAzimuth, currentRAAN[0], s=5, c='blue',                                 \
                     edgecolors='face' )

        ax5.scatter( currentAzimuth, currentAOP[0], s=5, c='blue',                                  \
                     edgecolors='face' )

        # ax6.scatter( e_cos_w[0], e_sin_w[0], s=5, c='blue',                                         \
        #              edgecolors='face' )

        # ax7.scatter( i_cos_raan[0], i_sin_raan[0], s=5, c='blue',                                   \
        #              edgecolors='face' )

# plotCircle( 1.0, ax6 )

ax1.set_xticks( np.arange( 0.0, 360.0, 30.0 ) )
ax1.set_xlim( 0.0, 360.0 )
ax1.grid( True )
ax1.set_xlabel( 'Launch azimuth [deg]' )
ax1.set_ylabel( 'Initial semi-major axis [m]' )
ax1.ticklabel_format( style='sci', axis='y', scilimits=(0,0), useoffset=False )
ax1.legend( markerscale=5 ).draggable( )

ax2.set_xticks( np.arange( 0.0, 360.0, 30.0 ) )
ax2.set_xlim( 0.0, 360.0 )
ax2.grid( True )
ax2.set_xlabel( 'Launch azimuth [deg]' )
ax2.set_ylabel( 'Initial eccentricity' )

ax3.set_xticks( np.arange( 0.0, 360.0, 30.0 ) )
ax3.set_xlim( 0.0, 360.0 )
ax3.grid( True )
ax3.set_xlabel( 'Launch azimuth [deg]' )
ax3.set_ylabel( 'Inclination [deg]' )

ax4.set_xticks( np.arange( 0.0, 360.0, 30.0 ) )
ax4.set_xlim( 0.0, 360.0 )
ax4.grid( True )
ax4.set_xlabel( 'Launch azimuth [deg]' )
ax4.set_ylabel( 'Right ascension of ascending node [deg]' )

ax5.set_xticks( np.arange( 0.0, 360.0, 30.0 ) )
ax5.set_xlim( 0.0, 360.0 )
ax5.grid( True )
ax5.set_xlabel( 'Launch azimuth [deg]' )
ax5.set_ylabel( 'Argument of periapsis [deg]' )

# ax6.set_ylabel('$e.sin(\omega)$')
# ax6.set_xlabel('$e.cos(\omega)$')
# ax6.grid(True)

# ax7.set_ylabel('$i.sin(\Omega)$')
# ax7.set_xlabel('$i.cos(\Omega)$')
# ax7.grid(True)

# Set global plot title
plt.suptitle( 'Initial orbital elements variation for launch velocity = '
               + str( current_velocity[ 0 ] ) + ' [m/s]' )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

# Show the plot
plt.show( )

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
