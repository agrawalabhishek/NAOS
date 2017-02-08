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

## Operations
# Connect to SQLite database.
try:
        database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity_with_perturbations/phase_45/leadingEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

phaseAngle = 45

## get data for escape cases
data1 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 );",                                      \
                     database )

data1.columns = [ 'vel_mag',                                                                \
                  'launch_azimuth' ]

escape_velocity_magnitude              = data1[ 'vel_mag' ]
escape_azimuth                         = data1[ 'launch_azimuth' ]

## get data for re-impact cases
data2 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( crash_flag = 1 );",                                       \
                     database )

data2.columns = [ 'vel_mag',                                                                \
                  'launch_azimuth' ]

crash_velocity_magnitude              = data2[ 'vel_mag' ]
crash_azimuth                         = data2[ 'launch_azimuth' ]

## get data for temporary capture cases
data3 = pd.read_sql( "SELECT    initial_velocity_x,                                         \
                                initial_velocity_y,                                         \
                                initial_velocity_z,                                         \
                                initial_inertial_velocity_x,                                \
                                initial_inertial_velocity_y,                                \
                                initial_inertial_velocity_z,                                \
                                ROUND( launch_azimuth ),                                    \
                                total_energy,                                               \
                                eccentricity                                                \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( end_flag = 1 );",                                         \
                     database )

data3.columns = [ 'rotFrame_vx',                                                            \
                  'rotFrame_vy',                                                            \
                  'rotFrame_vz',                                                            \
                  'inertial_vx',                                                            \
                  'inertial_vy',                                                            \
                  'inertial_vz',                                                            \
                  'launch_azimuth',                                                         \
                  'energy',                                                                 \
                  'eccentricity' ]

capture_rotFrame_vx                     = data3[ 'rotFrame_vx' ]
capture_rotFrame_vy                     = data3[ 'rotFrame_vy' ]
capture_rotFrame_vz                     = data3[ 'rotFrame_vz' ]
capture_inertial_vx                     = data3[ 'inertial_vx' ]
capture_inertial_vy                     = data3[ 'inertial_vy' ]
capture_inertial_vz                     = data3[ 'inertial_vz' ]
capture_azimuth                         = data3[ 'launch_azimuth' ]
capture_energy                          = data3[ 'energy' ]
capture_eccentricity                    = data3[ 'eccentricity' ]

if capture_azimuth.size != 0:
    ## Set up the figure
    fig = plt.figure( )
    gs  = gridspec.GridSpec( 3, 1, height_ratios = [ 1, 1, 1 ] )
    ax1 = plt.subplot( gs[ 0 ] )
    ax2 = plt.subplot( gs[ 1 ] )
    ax3 = plt.subplot( gs[ 2 ] )

    ## plot the histogram for number of particles escaped versus the initial launch velocity
    min_escape_velocity_magnitude = min(escape_velocity_magnitude)
    max_escape_velocity_magnitude = max(escape_velocity_magnitude)
    numberOfBins = (max_escape_velocity_magnitude - min_escape_velocity_magnitude)+1
    ax1.hist( escape_velocity_magnitude,
              bins=range(int(min_escape_velocity_magnitude), int(max_escape_velocity_magnitude)+2),
              align='left', color='green', alpha=0.8 )
    ax1.set_xlim( min_escape_velocity_magnitude, max_escape_velocity_magnitude )
    ax1.grid( )
else:
    ## Set up the figure
    # fig = plt.figure( )
    # gs = gridspec.GridSpec( 2, 1, height_ratios = [ 1, 1 ] )
    # ax1 = plt.subplot( gs[ 0 ] )
    # ax2 = plt.subplot( gs[ 1 ] )

    # ## plot the histogram for number of particles escaped versus the initial launch velocity
    # min_escape_velocity_magnitude = min(escape_velocity_magnitude)
    # max_escape_velocity_magnitude = max(escape_velocity_magnitude)
    # ax1.hist( escape_velocity_magnitude,
    #           bins=range(int(min_escape_velocity_magnitude), int(max_escape_velocity_magnitude)+2),
    #           align='left', color='green', alpha=0.8 )
    # ax1.set_xlim( min_escape_velocity_magnitude-1, max_escape_velocity_magnitude+1 )
    # ax1.grid( )

    # ## plot the histogram for number of particles that reimpacted versus initial launch velocity mag
    # min_crash_velocity_magnitude = min(crash_velocity_magnitude)
    # max_crash_velocity_magnitude = max(crash_velocity_magnitude)
    # ax2.hist( crash_velocity_magnitude,
    #           bins=range(int(min_crash_velocity_magnitude), int(max_crash_velocity_magnitude)+2),
    #           align='left', color='green', alpha=0.8 )
    # ax2.set_xlim( min_crash_velocity_magnitude, max_crash_velocity_magnitude )
    # ax2.grid( )

    ## plot the histogram for number of particles escaped and reimpacted versus
    ## initial launch velocity
    fig = plt.figure( )
    ax3 = fig.add_subplot( 111 )
    ax3.hist( [crash_velocity_magnitude, escape_velocity_magnitude],
              bins=range(1, 20+2),
              histtype='barstacked',
              color=['red', 'green'],
              label=['reimpact', 'escape'],
              align='left', alpha=0.7 )
    ax3.set_xlim( 0, 21 )
    ax3.xaxis.set_ticks( np.arange(0, 21, 1) )
    ax3.set_xlabel( '$V_{initial}$ [m/s]' )
    ax3.set_ylabel( 'Number of particles' )
    ax3.set_title( 'Particle fate for different initial launch velocities \n Phase angle = ' + str(phaseAngle) + '[deg]' )
    ax3.legend( )
    ax3.grid( )

## Show the plot
# plt.legend( handles = [ ax2Handle1, ax2Handle2, escapeBoundHandle, crashBoundHandle ],
#             bbox_to_anchor = ( 0.95, 0.5 ),
#             loc = 1 )
# plt.legend( bbox_to_anchor=(1.0, 1), loc=1, borderaxespad=0. )
plt.tight_layout( )
# plt.grid( )
plt.show( )
# plt.savefig('../data/trajectory_for_different_launch_azimuth/trajectoryPlot_'+str(figureIndex+1)+'.png', dpi=300)

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
