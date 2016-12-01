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
data = pd.read_csv( "../data/solutionParticleAroundUniformlyRotatingEllipsoid_orbitalElements.csv" )
cartesianData = pd.read_csv( "../data/solutionParticleAroundUniformlyRotatingEllipsoid.csv" )

semiMajor = data['semiMajor'].values
eccentricity = data['eccentricity'].values
inclination = data['inclination'].values
raan = data['raan'].values
aop = data['aop'].values
ta = data['ta'].values
t = data['time'].values

# convert time in seconds to earth days
t = t / ( 24.0 * 60.0 * 60.0 )

## Set up the figure
fig = plt.figure( )
ax1 = fig.add_subplot( 421 )
ax2 = fig.add_subplot( 422 )
ax3 = fig.add_subplot( 423 )
ax4 = fig.add_subplot( 424 )
ax5 = fig.add_subplot( 425 )

## plot the orbital elements
ax1.plot( t, semiMajor, color=colors.cnames['purple'] )
ax1.set_xlabel('time [Earth days]')
ax1.set_ylabel('semi-major axis [m]')
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
ax1.grid( )

ax2.plot( t, eccentricity, color=colors.cnames['purple'] )
ax2.set_xlabel('time [Earth days]')
ax2.set_ylabel('eccentricity')
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
ax2.grid( )

ax3.plot( t, inclination, color=colors.cnames['purple'] )
ax3.set_xlabel('time [Earth days]')
ax3.set_ylabel('inclination [deg]')
ax3.ticklabel_format(style='plain', useOffset=False)
ax3.grid( )

ax4.plot( t, raan, color=colors.cnames['purple'] )
ax4.set_xlabel('time [Earth days]')
ax4.set_ylabel('RAAN [deg]')
ax4.ticklabel_format(style='plain', useOffset=False)
ax4.grid( )

ax5.plot( t, aop, color=colors.cnames['purple'] )
ax5.set_xlabel('time [Earth days]')
ax5.set_ylabel('AOP [deg]')
ax5.ticklabel_format(style='plain', useOffset=False)
ax5.grid( )

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
