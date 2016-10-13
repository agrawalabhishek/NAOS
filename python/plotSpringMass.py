'''
Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
Distributed under the MIT License.
See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
'''

# Set up modules and packages
# I/O
import csv
from pprint import pprint

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
# Read data in csv file. data returned as a panda series.
data = pd.read_csv( '../data/springMassIntegrationSolution.csv' )
integrated_distance = data[ 'integrated_distance' ].values
integrated_velocity = data[ 'integrated_velocity' ].values
analytical_distance = data[ 'analytical_distance' ].values
analytical_velocity = data[ 'analytical_velocity' ].values
t = data[ 't' ].values

## Set up the figure
fig = plt.figure( )
ax1 = fig.add_subplot( 211 )
ax2 = fig.add_subplot( 212 )

## plot the distance
ax1.plot( t, integrated_distance, color=colors.cnames['red'], label='RK4 integrated' )
ax1.plot( t, analytical_distance, color=colors.cnames['green'], label='analytical' )
ax1.set_xlabel('time [s]')
ax1.set_ylabel('distance [m]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax1.grid( )
ax1.legend( )

## plot the velocity
ax2.plot( t, integrated_velocity, color=colors.cnames['red'], label='RK4 integrated' )
ax2.plot( t, analytical_velocity, color=colors.cnames['green'], label='analytical' )
ax2.set_xlabel('time [s]')
ax2.set_ylabel('velocity [m/s]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax2.grid( )
ax2.legend( )

## show plot
plt.tight_layout( )
plt.show( )

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
