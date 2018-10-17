
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D 

root_folder = "/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/"
folder = root_folder + "Examples/Ex63_3D_veins/"
input_folder = folder + "input_data/"
output_folder = folder + "output_data/"

#======================================================================
# FORWARD PROBLEM
#======================================================================
# Original
file_forward_signal = input_folder + "forwardSignal_98sensors.dat"
forward_signal = np.loadtxt(file_forward_signal)
# Code today
file_forward_signal_t = output_folder + "ForwardSignal.dat"
forward_signal_t = np.loadtxt(file_forward_signal)
# Difference
print(np.max(forward_signal-forward_signal_t))

#===========================
# TRAJECTORIES
#===========================
file_trajectories = output_folder + "Trajectory0.dat"
trajectories = np.loadtxt(file_trajectories)
nSteps, nRays = trajectories.shape
xCoord = trajectories[:, 1:nRays:3];
yCoord = trajectories[:, 2:nRays:3];
zCoord = trajectories[:, 3:nRays:3];
nRays = (int) nRays/3
# Figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
vecRays = np.arange(0, 3000, 1)
colors = cm.winter(np.linspace(0, 1, len(vecRays)))

for n, c in zip(vecRays, colors):
    plt.plot(xCoord[:, n], yCoord[:, n], zCoord[:, n], color=c)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
#ax.axis(np.array([0, 128*0.0001, 0, 128*0.0001, 0, 128*0.0001]))
plt.show(block=False)


#======================================================================
# ADJOINT PROBLEM
#======================================================================
# Original
file_adjoint_pressure = output_folder + "PixelPressure_original.dat"
adjoint_pressure = np.loadtxt(file_adjoint_pressure)

# Code today
file_adjoint_pressure_t = output_folder + "PixelPressure.dat"
adjoint_pressure_t = np.loadtxt(file_adjoint_pressure_t)

# Difference
print(np.mean(abs(adjoint_pressure-adjoint_pressure_t)))
print(np.mean(abs(adjoint_pressure)))
print(np.mean(abs(adjoint_pressure_t)))

for i in range(4):
    plt.figure()
    plt.imshow(adjoint_pressure[1+128*i:128*(i+1), :] - adjoint_pressure_t[1+128*i:128*(i+1), :])
    plt.colorbar()
    plt.show(block=False)
    
    plt.figure()
    plt.imshow(adjoint_pressure_t[1+128*i:128*(i+1), :])
    plt.colorbar()
    plt.show(block=False)
    
    plt.figure()
    plt.imshow(adjoint_pressure[1+128*i:128*(i+1), :])
    plt.colorbar()
    plt.show(block=False)
    
    input()
    plt.close()
    plt.close()
    plt.close()



