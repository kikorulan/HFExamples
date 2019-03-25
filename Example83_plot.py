
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
# PLOT
#================================================================================
#==============================
# Read data
#==============================
# Initial pressure
with h5py.File(data_folder + "photo_data_smooth.h5", "r") as f:
    photo_data = f["/Points"].value

photo_data = np.transpose(photo_data).copy()
# Load forward data - kWave
with h5py.File(data_folder + "forwardData_nImages_kWave.h5", "r") as f:
    forward_data_nImages_kWave = f["/forward_data"].value

#==============================
# Plot
#==============================
for index_volume in range(250, 254):
    # Load forward data
    with h5py.File(output_folder + "forwardData_vol" + str(index_volume+1) + "_python.h5", "r") as f:
        forward_data = f["/forward_data"].value
    
    forward_data_kWave = np.transpose(forward_data_nImages_kWave[index_volume, :, :])
    
    # Sensor signal
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(forward_data[200, :])
    plt.subplot(2, 1, 2)
    plt.plot(forward_data_kWave[200, :], c='r')
    plt.show(block=False)
    
    # Sinogram
    plt.figure()
    plt.imshow(forward_data, aspect='auto')
    plt.show(block=False)
    plt.figure()
    plt.imshow(forward_data_kWave, aspect='auto')
    plt.show(block=False)
    # Wait and close figures
    input("Press Enter to continue...")
    plt.close('all')

#============================
# Plot maximum projection
#============================
'''
# Proj XY
projXY = np.max(initial_pressure, axis=2)
plt.figure()
plt.imshow(projXY)
plt.show(block=False)
# Proj XZ
projXZ = np.max(initial_pressure, axis=1)
plt.figure()
plt.imshow(projXZ)
plt.show(block=False)
# Proj YZ
projYZ = np.max(initial_pressure, axis=0)
plt.figure()
plt.imshow(projYZ)
plt.colorbar()
plt.show(block=False)
'''
