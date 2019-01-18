##  % Read data from files

#==================================================
# LIBRARIES
#==================================================
import numpy as np
import os
import sys
sys.path.append('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/HighFreq_2DRT/utils_python')
import utils as u
# Path
path = '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex74_3D_thinveins'
os.chdir(path)

# Norm distance function
def norm_distance(x, y):
    return np.sqrt(np.sum(np.multiply(x.flatten()-y.flatten(), x.flatten()-y.flatten())))

#==================================================
# Dimensions
#==================================================
# Import dimensions
dim = np.loadtxt('input_data/dimensions.dat');
Nx = dim[0, 0]
Ny = dim[0, 1]
Nz = dim[0, 2] 
dx = dim[1, 0]
dy = dim[1, 1]
dz = dim[1, 2]
