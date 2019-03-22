
#============================================================
# FOLDER and FILE NAMES
#============================================================
# LIBRARIES
import numpy as np
import os
import socket
import sys
from datetime import datetime
import matplotlib.pyplot as plt

# EXAMPLE FOLDER
example = "Examples/Ex79_3D_norm/"
host_name = socket.gethostname()
if (host_name == "maryam.cs.ucl.ac.uk" or host_name == "ember.cs.ucl.ac.uk"):
    root_folder = "/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/"
    sys.path.append(root_folder + "HighFreq_3DRT/Build/bin")
elif (host_name == "hannover"):
    root_folder = "/home/wontek/sharedWK/Examples/"
    sys.path.append(root_folder + "HighFreq_3DRT/Build/bin")
else:
    root_folder = "/home/wonhong/sharedWK/"
    sys.path.append(root_folder + "RTlib")

example_folder = root_folder + example
# FOLDERS
input_folder = example_folder + "input_data/"
output_folder = example_folder + "output_data/"
os.chdir(example_folder)

# FILE NAMES
file_dimensions       = input_folder + "dimensions.dat"
file_initial_pressure = input_folder + "initial_pressure.dat"
file_sound_speed      = input_folder + "sound_speed.dat"
file_sensors_1a       = input_folder + "sensors_1a.dat"
file_sensors_1b       = input_folder + "sensors_1b.dat"
file_sensors_2        = input_folder + "sensors_2.dat"
file_sensors_10       = input_folder + "sensors_10.dat"
file_sensors_3600     = input_folder + "sensors_3600.dat"
file_sensors_14400    = input_folder + "sensors_14400.dat"

# SET GPU
import pyrtgrid as rt
GPU_INDEX = 1
rt.setGPU(GPU_INDEX)

#================================================================================
# FORWARD PROBLEM
#================================================================================
#==============================
# Load data
#==============================
# Dimensions
dim = np.loadtxt(file_dimensions)
# Sound speed
sound_speed = np.loadtxt(file_sound_speed)
sound_speed = np.reshape(sound_speed, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
sound_speed = np.transpose(sound_speed, (1, 2, 0)).copy()
# Initial pressure
initial_pressure = np.loadtxt(file_initial_pressure)
initial_pressure = np.reshape(initial_pressure, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
initial_pressure = np.transpose(initial_pressure, (1, 2, 0)).copy()
# Sensors
sensors = np.loadtxt(file_sensors_3600)
# Create object
rtobj = rt.rtgrid(dim)
# Load sound speed
rtobj.load_c(sound_speed)
# Load initial pressure
rtobj.load_u(initial_pressure)
# Load sensors
rtobj.load_sensors(sensors)
nSensors = 3600
sensor_array = np.arange(nSensors).astype(np.int64)

#==============================
# Plot initial pressure
#==============================
#projX = np.max(initial_pressure, axis=0)
#plt.figure()
#plt.imshow(projX)
#plt.colorbar()
#plt.show(block=False)

#==============================
# Apply operator
#==============================
factor_norm = 1e18
nIter = 10
norm_vector = np.zeros([nIter, 1])
adjoint_pressure_next = initial_pressure/factor_norm

# Iteration 
plt.close('all')
for i in range(nIter):
    print(i)
    adjoint_pressure_prev = np.array(factor_norm*adjoint_pressure_next, dtype=np.float64).copy()
    rtobj.load_u(adjoint_pressure_prev)
    forward_signal        = rtobj.forward_operator(sensor_array)
    adjoint_pressure_next = rtobj.adjoint_operator(sensor_array, forward_signal)
    norm_vector[i] = np.max(adjoint_pressure_next)/np.max(adjoint_pressure_prev)
    
    if i > nIter-3:
        plt.figure()
        plt.plot(forward_signal.transpose())
        plt.show(block=False)
        
        projY = np.max(adjoint_pressure_next, axis=1)
        plt.figure()
        plt.imshow(projY)
        plt.colorbar()
        plt.show(block=False)
        
        sliceY = adjoint_pressure_next[:, 81, :]
        plt.figure()
        plt.imshow(sliceY)
        plt.colorbar()
        plt.show(block=False)
 

plt.figure()
plt.plot(norm_vector)
plt.show(block=False)

plt.figure()
plt.plot(forward_signal.transpose())
plt.show(block=False)

projX = np.max(adjoint_pressure_next, axis=0)
plt.figure()
plt.imshow(projX)
plt.colorbar()
plt.show(block=False)


#==============================
# RESULTS
#==============================
tail_norm_vector = np.mean(norm_vector[20:])
print(tail_norm_vector)
# NORM 1 SENSOR  = 9.4e-19
# NORM 2 SENSORS = 9.4e-19
