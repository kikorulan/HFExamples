#==================================================================================================================================
# EXAMPLE 65 - WRAP IN PYTHON
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
folder = root_folder + "Examples/Ex65_python/"
input_folder = folder + "input_data/"
output_folder = folder + "output_data/"
os.chdir(folder)
nSensors = 98
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
sensor_array = np.arange(nSensors)
startTime = datetime.now() # time
forward_signal = rtobj.forward_operator(sensor_array)
endTime = datetime.now()
print(endTime - startTime)
 
#==============================
# Plot - forward
#==============================
# Comparison 
file_forward_signal = input_folder + "forwardSignal_98sensors.dat"
forward_signal_cpp = np.loadtxt(file_forward_signal)

for i in range(nSensors):
    #plt.figure()
    error = forward_signal[i, :] - forward_signal_cpp[i+1, :]
    plt.plot(error)

#plt.figure()
plt.show(block=False)


for i in range(5):
    plt.figure()
    plt.plot(forward_signal_cpp[i+1, :])
    plt.plot(forward_signal[i, :])
    plt.show(block=False)

#==============================
# Plot sound speed
#==============================
plt.figure()
plt.imshow(sound_speed[:, :, 1])
plt.colorbar()
plt.show(block=False)

#==============================
# Get filter
#==============================
##  filters = rtobj.get_filter()
##  plt.imshow(filters)
##  plt.show(block=False)

#================================================================================
# LOAD SIGNALS
#================================================================================
#==============================
# Get signal
#==============================
##  nS = 1
##  signal = np.random.random(1501).astype(np.float32)
##  signal_host = rtobj.get_signal(nS, signal)
##  plt.figure()
##  plt.plot(signal_host - signal)
##  plt.show(block=False)
##  plt.figure()
##  plt.plot(signal_host)
##  plt.show(block=False)
##  plt.figure()
##  plt.plot(signal)
##  plt.show(block=False)

#==============================
# Get signal convolve
#==============================
##  file_forward_signal = input_folder + "forwardSignal_98sensors.dat"
##  forward_signal = np.loadtxt(file_forward_signal)
##  nS = 1
##  signal = forward_signal[nS, :].astype(np.float32)
##  signal_convolve = rtobj.get_signalConvolve(nS, signal)
##  plt.figure()
##  plt.imshow(signal_convolve)
##  plt.show(block=False)

#================================================================================
# ADJOINT PROBLEM
#================================================================================
#==============================
# LOAD DATA: ADJOINT
#==============================
# Adjoint pressure
file_adjoint_pressure = input_folder + "adjoint_pressure.dat"
adjoint_pressure = np.loadtxt(file_adjoint_pressure)
adjoint_pressure = np.reshape(adjoint_pressure, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
adjoint_pressure = np.transpose(adjoint_pressure, (1, 2, 0)).copy()

#==============================
# Adjoint operator
#==============================
file_forward_signal = input_folder + "forwardSignal_98sensors.dat"
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

for i in range(4):
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

#===========================
# TRAJECTORIES
#===========================
##  file_trajectories_cpp = input_folder + "Trajectory0.dat"
##  trajectories_cpp = np.loadtxt(file_trajectories_cpp)
##  # Python
##  file_trajectories = output_folder + "Trajectory0.dat"
##  trajectories = np.loadtxt(file_trajectories)
##  # Difference
##  trajectories_nonan = trajectories
##  trajectories_nonan[np.isnan(trajectories)] = 0
##  trajectories_nonan[np.isinf(trajectories)] = 0
##  trajectories_cpp_nonan = trajectories_cpp
##  trajectories_cpp_nonan[np.isnan(trajectories_cpp)] = 0
##  trajectories_cpp_nonan[np.isinf(trajectories_cpp)] = 0
##  print(np.max(trajectories_nonan - trajectories_cpp_nonan))
##  # Trajectories
##  nSteps, nRays = trajectories.shape
##  xCoord = trajectories[:, 1:nRays:3]
##  yCoord = trajectories[:, 2:nRays:3]
##  zCoord = trajectories[:, 3:nRays:3]
##  nRays = nRays/3
##  # Figure
##  fig = plt.figure()
##  ax = fig.add_subplot(111, projection='3d')
##  vecRays = np.arange(0, 3000, 1)
##  colors = cm.winter(np.linspace(0, 1, len(vecRays)))
##  
##  for n, c in zip(vecRays, colors):
##      plt.plot(xCoord[:, n], yCoord[:, n], zCoord[:, n], color=c)
##  
##  ax.set_xlabel('X Label')
##  ax.set_ylabel('Y Label')
##  ax.set_zlabel('Z Label')
##  #ax.axis(np.array([0, 128*0.0001, 0, 128*0.0001, 0, 128*0.0001]))
##  plt.show(block=False)


#================================================================================
# ADJOINT SINGLE SENSOR 
#================================================================================
##  ref_folder = root_folder + "Examples/Ex63_3D_veins/output_data/"

#==============================
# LOAD DATA SINGLE SENSOR
#==============================
##  # Adjoint pressure
##  file_adjoint_sensor = ref_folder + "PixelPressure_" + str(nS) + ".dat"
##  adjoint_sensor = np.loadtxt(file_adjoint_sensor)
##  adjoint_sensor = np.reshape(adjoint_sensor, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
##  adjoint_sensor = np.transpose(adjoint_sensor, (1, 2, 0)).copy()

#==============================
# LOAD DATA SINGLE SENSOR
#==============================
##  # Adjoint pressure
##  file_adjoint_sensor_py = output_folder + "PixelPressure_" + str(nS) + ".dat"
##  adjoint_sensor = np.loadtxt(file_adjoint_sensor)
##  adjoint_sensor = np.reshape(adjoint_sensor, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
##  adjoint_sensor = np.transpose(adjoint_sensor, (1, 2, 0)).copy()


#==============================
# ADJOINT SINGLE SENSOR 0
#==============================
##  file_forward_signal = input_folder + "forwardSignal_98sensors.dat"
##  forward_signal = np.loadtxt(file_forward_signal).astype(np.float32)
##  adjoint_pressure_py = np.zeros((128, 128, 128)).astype(np.float32)
##  for nS in range(98):
##      forward_signal_sensor = forward_signal[nS+1:nS+2, :].copy()
##      sensor_array = np.array([nS]) 
##      startTime = datetime.now()
##      adjoint_sensor_py = rtobj.adjoint_operator(sensor_array, forward_signal_sensor)
##      adjoint_pressure_py = adjoint_pressure_py + adjoint_sensor_py
##      endTime = datetime.now()
##      print(endTime - startTime)
##  
##  
##  plt.figure()
##  plt.imshow(adjoint_pressure_py[:, :, 0])
##  plt.colorbar()
##  plt.show(block=False)
#==============================
# Plot
#==============================
##  print(np.mean(adjoint_sensor_py))
##  print(np.mean(adjoint_sensor))
##  print(np.mean(adjoint_sensor_py - adjoint_sensor))
##  
##  plt.figure()
##  plt.imshow(adjoint_sensor[:, :, i] - adjoint_sensor_py[:, :, i])
##  plt.colorbar()
##  plt.show(block=False)
##  
##  plt.figure()
##  plt.imshow(adjoint_sensor_py[:, :, i])
##  plt.colorbar()
##  plt.show(block=False)
##  
##  plt.figure()
##  plt.imshow(adjoint_sensor[:, :, i])
##  plt.colorbar()
##  plt.show(block=False)
##  
##  input()
##  plt.close()
##  plt.close()
##  plt.close()


