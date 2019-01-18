#==================================================================================================================================
# EXAMPLE 68 - WRAP IN PYTHON
#==================================================================================================================================
# Choose machine
import socket
host = socket.gethostname()
if (host == "hannover"):
    root_folder = "/home/wontek/sharedWK/"
elif (host == "ember.cs.ucl.ac.uk" or host == "maryam.cs.ucl.ac.uk"):
    root_folder = "/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/"
else:
    raise ValueError("Unknown machine")

# Add path
import sys
sys.path.append(root_folder + "HighFreq_3DRT/Build/bin")

# Import libraries
import numpy as np
import pyrtgrid as rt
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D 
import os

# Folder
folder = root_folder + "Examples/Ex68_3D_veins_resize/"
input_folder = folder + "input_data/"
output_folder = folder + "output_data/"
os.chdir(folder)
nSensors = 81
#================================================================================
# FORWARD PROBLEM
#================================================================================
#==============================
# Load data
#==============================
# Dimensions
file_dimensions = input_folder + "dimensions.dat"
dim = np.loadtxt(file_dimensions)
# Sound speed
file_sound_speed = input_folder + "sound_speed.dat"
sound_speed = np.loadtxt(file_sound_speed)
sound_speed = np.reshape(sound_speed, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
sound_speed = np.transpose(sound_speed, (1, 2, 0)).copy()
# Initial pressure
file_initial_pressure = input_folder + "initial_pressure_veins_80x240x240.dat"
initial_pressure = np.loadtxt(file_initial_pressure)
initial_pressure = np.reshape(initial_pressure, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
initial_pressure = np.transpose(initial_pressure, (1, 2, 0)).copy()
# Sensors
file_sensors = input_folder + "sensors_subsampled_81.dat"
sensors = np.loadtxt(file_sensors)
# Create object
rtobj = rt.rtgrid(dim)
# Load sound speed
rtobj.load_c(sound_speed)
# Load initial pressure
rtobj.load_u(initial_pressure)
# Load sensors
rtobj.load_sensors(sensors)

#==============================
# Forward operator
#==============================
sensor_array = np.arange(nSensors)
startTime = datetime.now() # time
forward_signal = rtobj.forward_operator(sensor_array)
endTime = datetime.now()
print(endTime - startTime)
 
#==============================
# Plot - forward
#==============================
# Comparison 
file_forward_signal = input_folder + "forwardSignal_reference_81sensors.dat"
forward_signal_cpp = np.loadtxt(file_forward_signal)

# Forward signal
plt.figure()
plt.imshow(forward_signal)
plt.colorbar()
plt.show(block=False)
# Error
errorImage = forward_signal - forward_signal_cpp[1:, :]
plt.figure()
plt.imshow(errorImage)
plt.colorbar()
plt.show(block=False)

# Errors for sensors
for i in range(5):
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(forward_signal_cpp[i+1, :])
    #plt.subplot(2, 1, 2)
    plt.plot(forward_signal[i, :], color='g')
    plt.subplot(2, 1, 2)
    plt.plot(forward_signal_cpp[i+1, :]-forward_signal[i, :], color='r')
    plt.legend()
    plt.show(block=False)

#================================================================================
# ADJOINT PROBLEM
#================================================================================
#==============================
# LOAD DATA: ADJOINT
#==============================
# Adjoint pressure
file_adjoint_pressure = input_folder + "pixelPressure_adjoint_81sensors.dat"
adjoint_pressure = np.loadtxt(file_adjoint_pressure)
adjoint_pressure = np.reshape(adjoint_pressure, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
adjoint_pressure = np.transpose(adjoint_pressure, (1, 2, 0)).copy()

#==============================
# Adjoint operator
#==============================
file_forward_signal = input_folder + "forwardSignal_reference_81sensors.dat"
forward_signal = np.loadtxt(file_forward_signal).astype(np.float32)
forward_signal = forward_signal[1:, :].copy()
sensor_array = np.arange(nSensors).astype(np.int32)
startTime = datetime.now() # time
adjoint_pressure_py = rtobj.adjoint_operator(sensor_array, forward_signal)
endTime = datetime.now()
print(endTime - startTime)

#==============================
# Plot
#==============================
print(np.max(adjoint_pressure_py))
print(np.max(adjoint_pressure))
print(np.max(adjoint_pressure_py - adjoint_pressure))

# Projections
projX_py = np.amax(adjoint_pressure_py, axis=0)
projY_py = np.amax(adjoint_pressure_py, axis=1)
projZ_py = np.amax(adjoint_pressure_py, axis=2)
projX = np.amax(adjoint_pressure, axis=0)
projY = np.amax(adjoint_pressure, axis=1)
projZ = np.amax(adjoint_pressure, axis=2)

# Python projections
plt.figure()
plt.imshow(projX_py)
plt.colorbar()
plt.show(block=False)

plt.figure()
plt.imshow(projY_py)
plt.colorbar()
plt.show(block=False)

plt.figure()
plt.imshow(projZ_py)
plt.colorbar()
plt.show(block=False)

# CPP projections
plt.figure()
plt.imshow(projX)
plt.colorbar()
plt.show(block=False)

plt.figure()
plt.imshow(projY)
plt.colorbar()
plt.show(block=False)

plt.figure()
plt.imshow(projZ)
plt.colorbar()
plt.show(block=False)
