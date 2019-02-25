
# LIBRARIES
import numpy as np
import socket
import os

#============================================================
# FOLDER and FILE NAMES
#============================================================
# EXAMPLE FOLDER
example = "Ex79_3D_norm/"
host_name = socket.gethostname()
if (host_name == "maryam.cs.ucl.ac.uk" or host_name == "ember.cs.ucl.ac.uk"):
    host_folder = "/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/"
elif (host_name == "hannover"):
    host_folder = "/home/wontek/sharedWK/Examples/"
else:
    host_folder = "/home/frullan/HighFreqCode/Examples/"
example_folder = host_folder + example
# FOLDERS
input_folder = example_folder + "input_data/"
output_folder = example_folder + "output_data/"
os.chdir(example_folder)

# FILE NAMES
file_name_dimensions       = input_folder + "dimensions.dat"
file_name_initial_pressure = input_folder + "initial_pressure.dat"
file_name_sound_speed      = input_folder + "sound_speed.dat"
file_name_sensors_1        = input_folder + "sensors_1.dat"
file_name_sensors_2        = input_folder + "sensors_2.dat"
file_name_sensors_14400    = input_folder + "sensors_14400.dat"

#============================================================
# GENERATE DOMAIN
#============================================================
Nx = 80
Ny = 240
Nz = 240
dx = 5.3e-5
dy = 5.3e-5
dz = 5.3e-5
# Write Dimensions
file_dimensions = open(file_name_dimensions, "w")
file_dimensions.write(str(Nx) + " " + str(Ny) + " " + str(Nz) + "\n" + str(dx) + " " + str(dy) + " " + str(dz))
file_dimensions.close()

#============================================================
# SOUND SPEED
#============================================================
c0 = 1582.01
sound_speed = c0*np.ones([Nx*Ny, Nz], dtype=np.float32)
np.savetxt(file_name_sound_speed, sound_speed, fmt='%.2f')

#============================================================
# INITIAL PRESSURE
#============================================================
np.random.seed(1)
lim_zero = 10
initial_pressure = np.zeros([Nx, Ny, Nz])
initial_pressure[lim_zero:-lim_zero, lim_zero:-lim_zero, lim_zero:-lim_zero] = np.random.rand(Nx-2*lim_zero, Ny-2*lim_zero, Nz-2*lim_zero)
initial_pressure = np.transpose(initial_pressure, (2, 0, 1)).copy()
initial_pressure = np.reshape(initial_pressure, (Nx*Nz, Ny)).copy()
np.savetxt(file_name_initial_pressure, initial_pressure, fmt='%.3f')

#============================================================
# SENSORS
#============================================================
nSensorsArray = 120
nRaysPhi = 1024 
nRaysTheta = 1024
dt = 1.667e-8
tMax = 3.32e-6

def str_prec(number, prec = '%.6f'):
    return str(prec % (number))

# ARRAY OF 14400 SENSORS
file_sensors_14400 = open(file_name_sensors_14400, "w")
file_sensors_14400.write(str(dt) + " " + str(tMax) + " 0 0 0 0 0 0 0\n")
for k in range(nSensorsArray):
    zPos = k*dz*(Nz-1)/(nSensorsArray-1)
    for i in range(nSensorsArray):
        yPos = i*dy*(Ny-1)/(nSensorsArray-1)
        file_sensors_14400.write("0 " + str_prec(yPos) + " " + str_prec(zPos) + " " + str(nRaysPhi) + " " + str(nRaysTheta) + " -1.57 1.57 0.04 3.1\n")
file_sensors_14400.close()

# 1 SENSOR
file_sensors_1 = open(file_name_sensors_1, "w")
file_sensors_1.write(str(dt) + " " + str(tMax) + " 0 0 0 0 0 0 0\n")
y1   = 40
z1   = 40
yPos = y1*dy*(Ny-1)/(nSensorsArray-1)
zPos = z1*dz*(Nz-1)/(nSensorsArray-1)
file_sensors_1.write("0 " + str_prec(yPos) + " " + str_prec(zPos) + " " + str(nRaysPhi) + " " + str(nRaysTheta) + " -1.57 1.57 0.04 3.1\n")
file_sensors_1.close()

# 2 SENSORS
file_sensors_2 = open(file_name_sensors_2, "w")
file_sensors_2.write(str(dt) + " " + str(tMax) + " 0 0 0 0 0 0 0\n")
y1   = 40
z1   = 40
yPos = y1*dy*(Ny-1)/(nSensorsArray-1)
zPos = z1*dz*(Nz-1)/(nSensorsArray-1)
file_sensors_2.write("0 " + str_prec(yPos) + " " + str_prec(zPos) + " " + str(nRaysPhi) + " " + str(nRaysTheta) + " -1.57 1.57 0.04 3.1\n")
y2   = 80
z2   = 60
yPos = y2*dy*(Ny-1)/(nSensorsArray-1)
zPos = z2*dz*(Nz-1)/(nSensorsArray-1)
file_sensors_2.write("0 " + str_prec(yPos) + " " + str_prec(zPos) + " " + str(nRaysPhi) + " " + str(nRaysTheta) + " -1.57 1.57 0.04 3.1\n")
file_sensors_2.close()
