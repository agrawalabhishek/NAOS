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

## Operations
# Connect to SQLite database.
try:
        # database = sqlite3.connect("../data/regolith_launched_from_leading_edge/multiple_launch_velocity/phase_0/simulation_time_9_months/leadingEdge.db")
        # database = sqlite3.connect( "../data/regolith_launched_from_longest_edge/spherical_asteroid/longestEdge.db" )
        database = sqlite3.connect("../data/regolith_launched_from_longest_edge/multiple_launch_velocity/simulation_time_9_months/longestEdge.db")

except sqlite3.Error, e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

velRange = tuple( range( 2, 17, 2 ) )

data1 = pd.read_sql( "SELECT    trajectory_id,                                                  \
                                position_x,                                                     \
                                position_y,                                                     \
                                position_z,                                                     \
                                ROUND( initial_velocity_magnitude ),                            \
                                ROUND( launch_azimuth )                                         \
                     FROM       regolith_trajectory_results                                     \
                     WHERE      ROUND( initial_velocity_magnitude ) IN " + str(velRange) + ";", \
                     database )

if database:
    database.close( )

data1.columns = [ 'trajectory_id',                                                          \
                  'x',                                                                      \
                  'y',                                                                      \
                  'z',                                                                      \
                  'init_vel_mag',                                                           \
                  'launch_azimuth' ]

trajectory_id       = data1[ 'trajectory_id' ]
x                   = data1[ 'x' ]
y                   = data1[ 'y' ]
z                   = data1[ 'z' ]
initial_velocity    = data1[ 'init_vel_mag' ]
azimuth             = data1[ 'launch_azimuth' ]

unique_trajectory_id = np.unique( trajectory_id )

fig = plt.figure( )
gs = gridspec.GridSpec( 1, 1 )
ax1 = plt.subplot( gs[ 0 ] )

print "Processing data now...\n"

mao_velocities      = []
hao_velocities      = []
uhao_velocities     = []
ehao_velocities     = []

for index in range( 0, len( unique_trajectory_id ) ):
    current_trajectory = unique_trajectory_id[ index ]
    current_data_indices = np.where( trajectory_id == current_trajectory )
    current_data_indices = current_data_indices[ 0 ]
    xTemp = x[ current_data_indices ]
    yTemp = y[ current_data_indices ]
    zTemp = z[ current_data_indices ]

    current_velocity = initial_velocity[ current_data_indices ]
    current_velocity = current_velocity.tolist( )
    current_velocity = current_velocity[ 0 ]

    rangeOfParticles = np.sqrt( (xTemp)**2 + (yTemp)**2 + (zTemp)**2 )

    # MAO
    MAO = np.where( (rangeOfParticles <= 3.0*alpha) & (rangeOfParticles > 2.0*alpha) )
    MAO = MAO[ 0 ]
    MAO = MAO.tolist()
    if MAO:
        mao_velocities.append(current_velocity)

    # HAO
    HAO = np.where( (rangeOfParticles <= 5.0*alpha) & (rangeOfParticles > 3.0*alpha) )
    HAO = HAO[ 0 ]
    HAO = HAO.tolist( )
    if HAO:
        hao_velocities.append(current_velocity)

    # UHAO
    UHAO = np.where( (rangeOfParticles <= 7.0*alpha) & (rangeOfParticles > 5.0*alpha) )
    UHAO = UHAO[ 0 ]
    UHAO = UHAO.tolist( )
    if UHAO:
        uhao_velocities.append(current_velocity)

    # EHAO
    EHAO = np.where( rangeOfParticles > 7.0*alpha )
    EHAO = EHAO[ 0 ]
    EHAO = EHAO.tolist( )
    if EHAO:
        ehao_velocities.append(current_velocity)

ax1.hist( [lao_mao_velocities, hao_velocities, uhao_velocities, ehao_velocities],       \
          bins=range(1, 16+2 ),                                                         \
          histtype='bar',                                                               \
          align='left',                                                                 \
          label=[ 'MAO', 'HAO', 'UHAO', 'EHAO' ] )
ax1.xaxis.set_ticks( np.arange(0, 16+1, 1) )
ax1.grid(True)
ax1.legend( ).draggable( )
ax1.set_ylabel( 'Regolith count' )
ax1.set_xlabel( 'Launch velocity [m/s]' )
ax1.set_title('Regolith count in different altitude regimes \n Ellipsoidal asteroid')

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
