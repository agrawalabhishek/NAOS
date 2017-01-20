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
# Read data in csv file. data returned as a panda series.
data = pd.read_csv( '../data/sun_asteroid_2BP/sunAsteroid2BP.csv' )

xBody           = data["x_body_frame"].values
yBody           = data["y_body_frame"].values
zBody           = data["z_body_frame"].values
vxBody          = data["vx_body_frame"].values
vyBody          = data["vy_body_frame"].values
vzBody          = data["vz_body_frame"].values
t               = data["t"].values
xInertial       = data["x_inertial_frame"].values
yInertial       = data["y_inertial_frame"].values
zInertial       = data["z_inertial_frame"].values
vxInertial      = data["vx_inertial_frame"].values
vyInertial      = data["vy_inertial_frame"].values
vzInertial      = data["vz_inertial_frame"].values
sma             = data["sma"].values
eccentricity    = data["eccentricity"].values
inclination     = data["inclination"].values
raan            = data["raan"].values
aop             = data["aop"].values
ta              = data["ta"].values
jacobian        = data["jacobi"].values
energy          = data["energy"].values

sunGravParameter = 1.32712440018 * 10.0e+20
regolithPosition = [ 25000.00, 0.0, 0.0 ]
sunPosition = [ xBody[200], yBody[200], zBody[200] ]
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
