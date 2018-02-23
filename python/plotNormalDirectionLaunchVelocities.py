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
    database = sqlite3.connect( "../data/guarantee_escape_speed/longest_edge/normalLaunch.db" )

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

## get the local directional escape speeds in rotating and inertial frame for each launch azimuth
data4 = pd.read_sql( "SELECT    directional_escape_speed,                                   \
                                directional_inertial_escape_speed,                          \
                                ROUND( initial_velocity_magnitude ),                        \
                                ROUND( launch_declination ),                                \
                                ROUND( launch_azimuth )                                     \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( escape_flag = 1 );",                                      \
                     database )

data4.columns = [ 'directional_escape_speed',                                               \
                  'directional_inertial_escape_speed',                                      \
                  'initial_velocity_magnitude',                                             \
                  'directional_escape_declination',                                         \
                  'directional_escape_azimuth' ]

directionalEscapeSpeed          = data4[ 'directional_escape_speed' ]
inertialDirectionalEscapeSpeed  = data4[ 'directional_inertial_escape_speed' ]
velMag                          = data4[ 'initial_velocity_magnitude' ]
directionalEscapeDeclination    = data4[ 'directional_escape_declination' ]
directionalEscapeAzimuth        = data4[ 'directional_escape_azimuth' ]

data5 = pd.read_sql( "SELECT    ROUND( initial_velocity_magnitude )                         \
                     FROM       regolith_trajectory_results                                 \
                     WHERE      ( crash_flag = 1 );",                                       \
                     database )

data5.columns = [ 'initial_velocity_magnitude' ]

crashVelMag = data5[ 'initial_velocity_magnitude' ]

if database:
    database.close( )

fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax4 = plt.subplot( gs[ 0 ] )

ax4.plot( velMag, label='Escape' )
ax4.plot( crashVelMag, label='Re-impact' )
ax4.plot( directionalEscapeSpeed, label='Conservative guaranted escape speed' )

ax4.grid( True )
ax4.set_ylabel( '$V_{Launch}$ [m/s]' )
ax4.legend( ).draggable( )
ax4.axes.xaxis.set_ticklabels([])

plt.show( )
