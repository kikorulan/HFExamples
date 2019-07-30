
# LIBRARIES
import numpy as np
import socket
import os
import getpass

#============================================================
# FOLDER and FILE NAMES
#============================================================
# EXAMPLE FOLDER
example = "Ex84_ROI/"
user_name = getpass.getuser()
host_name = socket.gethostname()
if (host_name == "maryam.cs.ucl.ac.uk" or host_name == "ember.cs.ucl.ac.uk"):
    host_folder = "/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/"
elif (host_name == "hannover"):
    host_folder = "/home/wontek/sharedWK/Examples/"
else:
    host_folder = "/home/wonhong/sharedWK/Examples/"
        
example_folder = host_folder + example
# FOLDERS
input_folder = example_folder + "input_data/"
output_folder = example_folder + "output_data/"
os.chdir(example_folder)

#============================================================
# CONFIGURATION
#============================================================
# DOMAIN
Nx = 80
Ny = 240
Nz = 240
dx = 5.3e-5
dy = 5.3e-5
dz = 5.3e-5
# SOUND SPEED
c0 = 1580.0
# ROI
#x_roi_rel = 0.5
#y_roi_rel = 0.83
#z_roi_rel = 0.5
rad_roi_pixels = np.ceil(1e-3/dx)
x_roi = 3.5e-3 # Nx*dx*x_roi_rel
y_roi = 3.5e-3 # Ny*dy*y_roi_rel
z_roi = 9.8e-3 # Nz*dz*z_roi_rel
rad_roi = rad_roi_pixels*dx
# Step and max time
dt   = 1.667e-8
tMax = 1.2e-05
# Number of Rays
nRaysPhi = 100
nRaysTheta = 100
# Sensor array
nSensorsArray = 6

#============================================================
# WRITE TO FILE
#============================================================
# File names
file_name_dimensions  = input_folder + "dimensions.dat"
file_name_sound_speed = input_folder + "sound_speed.dat"
file_name_sensors_ROI = input_folder + "sensors_ROI_" + str(nSensorsArray*nSensorsArray) + ".dat"
file_name_pressure0   = input_folder + "pixelPressure_0.dat"

# Dimensions
file_dimensions = open(file_name_dimensions, "w")
file_dimensions.write(str(Nx) + " " + str(Ny) + " " + str(Nz) + "\n" + str(dx) + " " + str(dy) + " " + str(dz))
file_dimensions.close()

# Sound Speed
sound_speed = c0*np.ones([Nx*Ny, Nz], dtype=np.float32)
np.savetxt(file_name_sound_speed, sound_speed, fmt='%.2f')

# Pixel 0
pressure0 = np.zeros([Nx*Ny, Nz], dtype=np.float32)
np.savetxt(file_name_pressure0, pressure0, fmt='%.2f')

# Region of Interest and Sensors
def str_prec(number, prec = '%.6f'):
    return str(prec % (number))
file_sensors_ROI = open(file_name_sensors_ROI, "w")
file_sensors_ROI.write(str_prec(x_roi) + " " + str_prec(y_roi) + " " + str_prec(z_roi) + " " + str_prec(rad_roi) + "\n")
file_sensors_ROI.write(str(dt) + " " + str(tMax) + " 0 0\n")
file_sensors_ROI.write(str(nSensorsArray*nSensorsArray) + " " + str(nRaysPhi) + " " + str(nRaysTheta) + " 0\n")
for k in range(nSensorsArray):
    zPos = k*dz*(Nz-1)/(nSensorsArray-1)
    for i in range(nSensorsArray):
        yPos = i*dy*(Ny-1)/(nSensorsArray-1)
        file_sensors_ROI.write("0 " + str_prec(yPos) + " " + str_prec(zPos) + " 0\n")

file_sensors_ROI.close()

