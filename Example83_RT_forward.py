
#==================================================================================================================================
# EXAMPLE 83 - WRAP IN PYTHON
#==================================================================================================================================

# Import Libraries
import socket
import sys
import numpy as np
from datetime import datetime
import os
import h5py
import matplotlib.pyplot as plt
# Choose machine
host = socket.gethostname()
if (host == "hannover"):
    root_folder = "/home/wontek/sharedWK/"
    sys.path.append(root_folder + "HighFreq_3DRT/Build/bin")
elif (host == "ember.cs.ucl.ac.uk" or host == "maryam.cs.ucl.ac.uk" or host == "blaze.cs.ucl.ac.uk"):
    root_folder = "/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/"
    data_folder = "/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/MLproject/data/"
    sys.path.append(root_folder + "HighFreq_3DRT/Build/bin")
else:
    root_folder = "/home/wonhong/sharedWK/"
    data_folder = "/home/wonhong/sharedWK/training/"
    sys.path.append(root_folder + "RTlib")

# Choose GPU
GPU_INDEX = 1

# Import libraries
import pyrtgrid as rt

# Folder
folder = root_folder + "Examples/Ex83_python_3600/"
input_folder = folder + "input_data/"
output_folder = folder + "output_data/"
os.chdir(folder)

#================================================================================
# FORWARD PROBLEM
#================================================================================
# Input files
file_dimensions  = input_folder + "dimensions.dat"
file_sound_speed = input_folder + "sound_speed.dat"
file_sensors     = input_folder + "sensors_3600.dat"
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
with h5py.File(data_folder + "photo_data_smooth.h5", "r") as f:
    photo_data = f["/Points"].value

photo_data = np.transpose(photo_data).copy()
# Sensors
sensors = np.loadtxt(file_sensors)
# SET GPU
rt.setGPU(GPU_INDEX)
# Create object
rtobj = rt.rtgrid(dim)
# Load sound speed
rtobj.load_c(sound_speed)
# Load sensors
nSensors = 3600
sensor_array = np.arange(nSensors).astype(np.int64)
rtobj.load_sensors(sensors)

#==================================================
# LOOP OVER VOLUMES
#==================================================
for index_volume in range(250, 500):
    # Load initial pressure
    initial_pressure = photo_data[index_volume, :, :, :]
    rtobj.load_u(initial_pressure)
    
    # Forward operator
    startTime = datetime.now()
    forward_data = rtobj.forward_operator(sensor_array)
    endTime = datetime.now()
    print(endTime - startTime)
    
    # Write to disk
    hf = h5py.File(output_folder + "forwardData_vol" + str(index_volume+1) + "_python.h5", "w")
    hf.create_dataset('forward_data', data=forward_data)
    hf.close()

#================================================================================
# PLOT
#================================================================================
#==============================
# Read data
#==============================
##  with h5py.File(output_folder + "forwardData_vol" + str(index_volume+1) + "_kWave.h5", "r") as f:
##      forward_data_kWave = f["/forward_data"].value
##  
##  plt.figure()
##  plt.subplot(2, 1, 1)
##  plt.plot(forward_data[0, :])
##  plt.subplot(2, 1, 2)
##  plt.plot(forward_data_kWave[:, 0], c='r')
##  plt.show(block=False)
##  
##  plt.figure()
##  plt.imshow(forward_data, aspect='auto')
##  plt.show(block=False)
##  plt.figure()
##  plt.imshow(np.transpose(forward_data_kWave), aspect='auto')
##  plt.show(block=False)

#============================
# Plot maximum projection
#============================
##  # Proj XY
##  projXY = np.max(initial_pressure, axis=2)
##  plt.figure()
##  plt.imshow(projXY)
##  plt.show(block=False)
##  # Proj XZ
##  projXZ = np.max(initial_pressure, axis=1)
##  plt.figure()
##  plt.imshow(projXZ)
##  plt.show(block=False)
##  # Proj YZ
##  projYZ = np.max(initial_pressure, axis=0)
##  plt.figure()
##  plt.imshow(projYZ)
##  plt.colorbar()
##  plt.show(block=False)
