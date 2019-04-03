
# LIBRARIES
import numpy as np
import socket
import os
import getpass

#============================================================
# FOLDER and FILE NAMES
#============================================================
# EXAMPLE FOLDER
example = "Ex83_python_3600/"
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

# FILE NAMES
file_name_dimensions       = input_folder + "dimensions.dat"
file_name_sound_speed      = input_folder + "sound_speed.dat"
file_name_sensors_3600     = input_folder + "sensors_3600.dat"

#============================================================
# GENERATE DOMAIN
#============================================================
Nx = 40
Ny = 120
Nz = 120
dx = 1e-4
dy = 1e-4
dz = 1e-4
# Write Dimensions
file_dimensions = open(file_name_dimensions, "w")
file_dimensions.write(str(Nx) + " " + str(Ny) + " " + str(Nz) + "\n" + str(dx) + " " + str(dy) + " " + str(dz))
file_dimensions.close()

#============================================================
# SOUND SPEED
#============================================================
c0 = 1500.0
sound_speed = c0*np.ones([Nx*Ny, Nz], dtype=np.float32)
np.savetxt(file_name_sound_speed, sound_speed, fmt='%.2f')

#============================================================
# SENSORS
#============================================================
nRaysPhi   = 512
nRaysTheta = 512
dt   = 3e-8
Nt   = 250
tMax = dt*(Nt-0.5)

def str_prec(number, prec = '%.6f'):
    return str(prec % (number))

# ARRAY OF 3600 SENSORS
nSensorsArray = 60
file_sensors_3600 = open(file_name_sensors_3600, "w")
file_sensors_3600.write(str(dt) + " " + str(tMax) + " 0 0 0 0 0 0 0\n")
for k in range(nSensorsArray):
    zPos = k*dz*(Nz-1)/(nSensorsArray-1)
    for i in range(nSensorsArray):
        yPos = i*dy*(Ny-1)/(nSensorsArray-1)
        file_sensors_3600.write("0 " + str_prec(yPos) + " " + str_prec(zPos) + " " + str(nRaysPhi) + " " + str(nRaysTheta) + " -1.57 1.57 0.04 3.1\n")

file_sensors_3600.close()

