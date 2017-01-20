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

sunGravParameter = 1.32712440018 * 10.0e+20
regolithPosition = [ 25000.00, 0.0, 0.0 ]
sunPosition = [ 1.60463e+11, -5.48052e+10, 3.67516e+08 ]
sunPositionMagnitude = np.linalg.norm( sunPosition )

sunToRegolithPosition = np.subtract( regolithPosition, sunPosition )
sunToRegolithPositionMagnitude = np.linalg.norm( sunToRegolithPosition )

thirdBodyAccelerationTermOne = np.multiply( ( -1.0 * sunGravParameter / sunToRegolithPositionMagnitude**3 ), sunToRegolithPosition )
thirdBodyAccelerationTermTwo = np.multiply( ( -1.0 * sunGravParameter / sunPositionMagnitude**3 ), sunPosition )

thirdBodyPerturbation = np.add( thirdBodyAccelerationTermOne, thirdBodyAccelerationTermTwo )

print "third body perturbation acceleration:"
print thirdBodyPerturbation
print " "

rho = 1.0
solarConstant = 1.0e17
regolithGrainDensity = 3.2 * 1.0e3
regolithGrainRadius = 40.0e-6
areaToMassRatio = 3.0 / ( 4.0 * regolithGrainRadius * regolithGrainDensity )
multiplicationConstant = -1.0 * ( 1.0 + rho ) * solarConstant * areaToMassRatio

regolithToSunPosition = np.subtract( sunPosition, regolithPosition )
regolithToSunPositionMagnitude = np.linalg.norm( regolithToSunPosition )

srpAcceleration = np.multiply( ( multiplicationConstant / regolithToSunPositionMagnitude**3 ), regolithToSunPosition )

print "SRP perturbation acceleration"
print srpAcceleration

# Stop timer
end_time = time.time( )

# Print elapsed time
print "Script time: " + str("{:,g}".format(end_time - start_time)) + "s"

print ""
print "------------------------------------------------------------------"
print "                         Exited successfully!                     "
print "------------------------------------------------------------------"
print ""
