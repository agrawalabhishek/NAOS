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
data = pd.read_csv( '../data/solutionParticleAroundRotatingSpheroid_orbitalElements.csv' )

t = data[ 'time' ].values
energy = data[ 'energy' ].values
angularMomentum = data[ 'angularMomentum' ].values

# convert time in seconds to earth days
t = t / ( 24.0 * 60.0 * 60.0 )

## Set up the figure
fig = plt.figure( )
ax1 = fig.add_subplot( 211 )
ax2 = fig.add_subplot( 212 )

## plot the energy
ax1.plot( t, energy, color=colors.cnames['purple'] )
ax1.set_xlabel('time [Earth days]')
ax1.set_ylabel('energy')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax1.grid( )

## plot angular momentum
ax2.plot( t, angularMomentum, color=colors.cnames['purple'] )
ax2.set_xlabel('time [Earth days]')
ax2.set_ylabel('angular momentum')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useOffset=False)
ax2.grid( )

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
