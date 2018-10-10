#==================================================================================================================================
# EXAMPLE 65 - WRAP IN PYTHON
#==================================================================================================================================
# Add path
import sys
sys.path.append("/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/HighFreq_3DRT/Build/bin")
# Import libraries
import numpy as np
import pyrtgrid as rt
from datetime import datetime
import matplotlib.pyplot as plt
# Folder
folder = "/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex65_python/"
input_folder = folder + "input_data/"
output_folder = folder + "output_data/"

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
file_initial_pressure = input_folder + "initial_pressure_veins.dat"
initial_pressure = np.loadtxt(file_initial_pressure)
initial_pressure = np.reshape(initial_pressure, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
initial_pressure = np.transpose(initial_pressure, (1, 2, 0)).copy()
# Sensors
file_sensors = input_folder + "sensors.dat"
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
sensor_array = np.arange(98)
startTime = datetime.now() # time
forward_signal = rtobj.forward_operator(sensor_array)
endTime = datetime.now()
print(endTime - startTime)
 
#==============================
# Plot
#==============================
# Comparison 
file_forward_signal = input_folder + "forwardSignal_98sensors.dat"
forward_signal_cpp = np.loadtxt(file_forward_signal)

for i in range(98):
    #plt.figure()
    error = forward_signal[i, :] - forward_signal_cpp[i+1, :]
    plt.plot(error)

#plt.figure()
#plt.show(block=False)


#================================================================================
# ADJOINT PROBLEM
#================================================================================
#==============================
# Load data
#==============================
# Adjoint pressure
file_adjoint_pressure = input_folder + "adjoint_pressure.dat"
adjoint_pressure = np.loadtxt(file_adjoint_pressure)
adjoint_pressure = np.reshape(adjoint_pressure, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
adjoint_pressure = np.transpose(adjoint_pressure, (1, 2, 0)).copy()

#==============================
# Adjoint operator
#==============================
sensor_array = np.arange(98)
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



for i in range(10):
    plt.figure()
    plt.imshow(adjoint_pressure[:, :, i] - adjoint_pressure_py[:, :, i])
    plt.colorbar()
    plt.show(block=False)
    
    plt.figure()
    plt.imshow(adjoint_pressure_py[:, :, i])
    plt.colorbar()
    plt.show(block=False)
    
    plt.figure()
    plt.imshow(adjoint_pressure[:, :, i])
    plt.colorbar()
    plt.show(block=False)
    
    input()
    plt.close()
    plt.close()
    plt.close()
