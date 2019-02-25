
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
# LOAD INFO
#============================================================
