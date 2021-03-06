#==================================================================================================================================
# EXAMPLE 65 - WRAP IN PYTHON
#==================================================================================================================================
# Choose machine
import socket
import sys
host = socket.gethostname()
if (host == "hannover"):
    root_folder = "/home/wontek/sharedWK/"
    sys.path.append(root_folder + "HighFreq_3DRT/Build/bin")
elif (host == "ember.cs.ucl.ac.uk" or host == "maryam.cs.ucl.ac.uk"):
    root_folder = "/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/"
    sys.path.append(root_folder + "HighFreq_3DRT/Build/bin")
else:
    root_folder = "/home/wonhong/sharedWK/"
    sys.path.append(root_folder + "RTlib")

# Choose GPU
GPU_INDEX = 0
# Import libraries
import numpy as np
#import pyrtgrid as rt
import pyrtgrid_dummy as rt_dummy
from datetime import datetime

import os

# Folder
folder = root_folder + "Examples/Ex76_3D_40x120x120/"
input_folder = folder + "input_data/"
output_folder = folder + "output_data/"
os.chdir(folder)
nSensors = 36

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
file_initial_pressure = input_folder + "initial_pressure_veins_40x120x120.dat"
initial_pressure = np.loadtxt(file_initial_pressure)
initial_pressure = np.reshape(initial_pressure, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
initial_pressure = np.transpose(initial_pressure, (1, 2, 0)).copy()
# Sensors
file_sensors = input_folder + "sensors_subsampled_36.dat"
sensors = np.loadtxt(file_sensors)
# SET GPU
rt_dummy.setGPU(GPU_INDEX)
# Create object
rtobj = rt_dummy.rtgrid(dim)
##  # Load sound speed
rtobj.load_c(sound_speed)
# Load initial pressure
rtobj.load_u(initial_pressure)
# Load sensors
rtobj.load_sensors(sensors)
  
#==============================
# Forward operator
#==============================
sensor_array = np.arange(nSensors).astype(np.int64)
startTime = datetime.now() # time
forward_signal = rtobj.forward_operator(sensor_array)
endTime = datetime.now()
print(endTime - startTime)
 
#==============================
# Plot - forward
#==============================
# Comparison 
file_forward_signal = input_folder + "forwardSignal_reference_36sensors.dat"
forward_signal_cpp = np.loadtxt(file_forward_signal)


print(np.max(forward_signal))
print(np.max(forward_signal_cpp[1:, :]))
print(np.max(forward_signal - forward_signal_cpp[1:, :]))


#================================================================================
# ADJOINT PROBLEM
#================================================================================
#==============================
# LOAD DATA: ADJOINT
#==============================
# Adjoint pressure
file_adjoint_pressure = input_folder + "adjoint_pressure.dat"
adjoint_pressure_cpp = np.loadtxt(file_adjoint_pressure)
adjoint_pressure_cpp = np.reshape(adjoint_pressure_cpp, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
adjoint_pressure_cpp = np.transpose(adjoint_pressure_cpp, (1, 2, 0)).copy()

#==============================
# Adjoint operator
#==============================
file_forward_signal = input_folder + "forwardSignal_reference_36sensors.dat"
forward_signal = np.loadtxt(file_forward_signal).astype(np.float32)
forward_signal = forward_signal[1:, :].copy()
sensor_array = np.arange(nSensors).astype(np.int64)
startTime = datetime.now() # time
adjoint_pressure = rtobj.adjoint_operator(sensor_array, forward_signal)
endTime = datetime.now()
print(endTime - startTime)

#==============================
# Plot
#==============================
print(np.max(adjoint_pressure))
print(np.max(adjoint_pressure_cpp))
print(np.max(adjoint_pressure - adjoint_pressure_cpp))

print(np.mean(adjoint_pressure))
print(np.mean(adjoint_pressure_cpp))
print(np.mean(adjoint_pressure - adjoint_pressure_cpp))
