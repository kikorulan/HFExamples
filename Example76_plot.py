
#============================================================
# Import Libraries
#============================================================
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
GPU_INDEX = 1
# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import os

# Folder
folder = root_folder + "Examples/Ex76_3D_40x120x120/"
input_folder = folder + "input_data/"
output_folder = folder + "output_data/"
os.chdir(folder)

#======================================================================
# Read forward data
#======================================================================
# Dimensions
file_dimensions = input_folder + "dimensions.dat"
dim = np.loadtxt(file_dimensions)

# Initial pressure
file_initial_pressure = input_folder + "initial_pressure_veins_40x120x120.dat"
initial_pressure = np.loadtxt(file_initial_pressure)
initial_pressure = np.reshape(initial_pressure, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
initial_pressure = np.transpose(initial_pressure, (1, 2, 0)).copy()

# Forward signal - reference
file_forwardSignal_reference = input_folder + "forwardSignal_reference_1600sensors.dat"
forwardSignal_reference = np.loadtxt(file_forwardSignal_reference)
forwardSignal_reference = forwardSignal_reference[2:, :]

# Forward signal 
file_forwardSignal = output_folder + "ForwardSignal.dat"
forwardSignal = np.loadtxt(file_forwardSignal)
forwardSignal = forwardSignal[2:, :]

# Difference
print(np.max(forwardSignal_reference))
print(np.max(forwardSignal))
print(np.max(forwardSignal-forwardSignal_reference))

#==============================
# ADJOINT
#==============================
# Adjoint pressure - reference
file_adjoint_reference = input_folder + "pixelPressure_adjoint_1600sensors.dat"
adjoint_pressure_reference = np.loadtxt(file_adjoint_reference)
adjoint_pressure_reference = np.reshape(adjoint_pressure_reference, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
adjoint_pressure_reference = np.transpose(adjoint_pressure_reference, (1, 2, 0)).copy()
# Adjoint pressure
file_adjoint = output_folder + "PixelPressure.dat"
adjoint_pressure = np.loadtxt(file_adjoint)
adjoint_pressure = np.reshape(adjoint_pressure, (int(dim[0, 2]), int(dim[0, 0]), int(dim[0, 1])))
adjoint_pressure = np.transpose(adjoint_pressure, (1, 2, 0)).copy()

# Difference
print(np.max(adjoint_pressure_reference))
print(np.max(adjoint_pressure))
print(np.max(adjoint_pressure-adjoint_pressure_reference))


# XY projection
XYproj_ref = np.max(adjoint_pressure_reference, axis=2)
plt.figure()
plt.imshow(XYproj_ref)
plt.show(block=False)
XYproj = np.max(adjoint_pressure, axis=2)
plt.figure()
plt.imshow(XYproj)
plt.show(block=False)
