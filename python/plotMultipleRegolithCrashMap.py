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
beginRange = 0
endRange = 360

alpha = 20000.0
beta = 7000.0
gamma = 7000.0

## Set up the figure
fig = plt.figure( )
ax1 = fig.add_subplot(111)

## data in body frame
for num in range( beginRange, endRange, 10 ):
    dynamicString = '../data/trajectory_for_different_launch_azimuth/'
    dynamicString = dynamicString + 'regolithTrajectoryAtAzimuth'
    dynamicString = dynamicString + str(num) + '.csv'

    data = pd.read_csv( dynamicString )
    x = data[ 'x' ].values
    y = data[ 'y' ].values
    z = data[ 'z' ].values
    vx = data[ 'vx' ].values
    vy = data[ 'vy' ].values
    vz = data[ 'vz' ].values
    t = data[ 't' ].values

    # check whether the final coordinate is on the surface of asteroid or not
    finalPoint = len( x ) - 1
    crashCheck = x[ finalPoint ]**2 / alpha**2       \
                    + y[ finalPoint ]**2 / beta**2   \
                    + z[ finalPoint ]**2 / gamma**2  \
                    - 1.0

    if abs( crashCheck ) <= 1.0e-12:
        c = np.random.random( )
        plotColor = cm.rainbow( num )

        ## calculate the lat long for starting point
        startRadialDistance = np.sqrt( x[0]**2 + y[0]**2 + z[0]**2 )
        startLongitude = np.arctan2( y[0], x[0] ) * 180.0 / np.pi
        startLatitude = np.arcsin( z[0] / startRadialDistance ) * 180 / np.pi

        ## calculate lat long for end points
        endRadialDistance = np.sqrt( x[finalPoint]**2 + y[finalPoint]**2 + z[finalPoint]**2 )
        endLongitude = np.arctan2( y[finalPoint], x[finalPoint] ) * 180 / np.pi
        endLatitude = np.arcsin( z[finalPoint] / endRadialDistance ) * 180 / np.pi

        ## indicate starting point
        ax1.scatter( startLongitude, startLatitude, color='black' )
        ax1.text( startLongitude, startLatitude, 'start', size=10, zorder=1, color='black' )

        ## Plot end points on lat long map
        labelString = 'Azimuth = ' + str( num )
        ax1.scatter( endLongitude, endLatitude, color=plotColor, label=labelString )

        ## indicate ending point
        ax1.text( endLongitude, endLatitude, 'end', size=10, zorder=1, color=plotColor )
    else:
        continue

## format axis and title
ax1.set_xlabel('longitude [degree]')
ax1.set_ylabel('latitude [degree]')
ax1.set_xlim( -180.0, 180.0 )
ax1.set_ylim( -90.0, 90.0 )
# ax1.ticklabel_format(style='plain', axis='both', useoffset=False)
ax1.set_title( 'Regolith crash map for asteroid Eros' )
ax1.legend( loc='center left', bbox_to_anchor=( 0.995,0.5 ) )

## Show the plot
# plt.legend( )
# plt.tight_layout( )
plt.grid( )
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
