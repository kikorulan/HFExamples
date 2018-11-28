#!/bin/bash

#================================================================================
# EXAMPLE 68
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================

# Output folder
export EXAMPLE_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex70_3D_synchronization/"
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"03/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
#export INITIAL_PRESSURE="initial_pressure_veins_80x240x240.dat"
export INITIAL_PRESSURE="initial_pressure_veins_smooth.dat"
export SENSORS="sensors_subsampled_1sensor.dat" 
export FORWARD_SIGNAL="forwardSignal_delay_1_0_delta.dat"

# Mode
export MODE="-a"

# Generate dimensions file
Nx=80  dx=0.000053
Ny=240 dy=0.000053
Nz=240 dz=0.000053
cat > $INPUT_FOLDER$DIMENSIONS <<EOF
$Nx $Ny $Nz
$dx $dy $dz
EOF

# Parameters
nSensorsArray=5
nRaysPhi=1024
nRaysTheta=1024
dt=1.6667e-8
tMax=8.0836e-06
##  # Generate sensor file
##  cat > $INPUT_FOLDER$SENSORS<<EOF
##  EOF
##  # Step and tMax
##  echo "$dt $tMax 0 0 0 0 0 0 0" >> $INPUT_FOLDER$SENSORS
##  # YZ
##  for ((k=0; k<nSensorsArray; k++)); do
##      zPos=$(echo "scale=4;($k*$dz*($Nz-1))/($nSensorsArray-1)" | bc)
##      for ((i=0; i<nSensorsArray; i++)); do
##          yPos=$(echo "scale=4;($i*$dy*($Ny-1))/($nSensorsArray-1)" | bc)
##          echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.04 3.1" >> $INPUT_FOLDER$SENSORS
##      done 
##  done 

# Generate sensor file
cat > $INPUT_FOLDER$SENSORS<<EOF
EOF
echo "$dt $tMax 0 0 0 0 0 0 0" >> $INPUT_FOLDER$SENSORS
yPos=$(echo "scale=6;($dy*$Ny/2)" | bc)
zPos=$(echo "scale=6;($dy*$Nz/2)" | bc)
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.04 3.1" >> $INPUT_FOLDER$SENSORS

# Call RT solver
export OMP_NUM_THREADS=26

#====================
# CUDA MEMCHECK
#====================
#/opt/cuda/cuda-8.0/bin/cuda-memcheck \
RTsolver_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                                                  $INPUT_FOLDER$INITIAL_PRESSURE $INPUT_FOLDER$SENSORS \
                                                  $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL
