#!/bin/bash

#================================================================================
# EXAMPLE 57
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================

# Output folder
export EXAMPLE_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex67_3D_veins_homo/"
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"03/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export INITIAL_PRESSURE="initial_pressure_veins.dat"
export SENSORS="sensors.dat" 
export FORWARD_SIGNAL="forwardSignal_98sensors.dat"

# Mode
export MODE="-a"

# Generate dimensions file
Nx=128 dx=0.0001
Ny=128 dy=0.0001
Nz=128 dz=0.0001
cat > $INPUT_FOLDER$DIMENSIONS <<EOF
$Nx $Ny $Nz
$dx $dy $dz
EOF

# Parameters
nSensorsArray=5
nRaysPhi=1000  # 1000 
nRaysTheta=300 # 300
dt=1e-8
tMax=1.5e-5
# Generate sensor file
cat > $INPUT_FOLDER$SENSORS<<EOF
EOF
# Step and tMax
echo "$dt $tMax 0 0 0 0 0 0 0" >> $INPUT_FOLDER$SENSORS
# XY Bottom
zPos=$(echo "scale=4;($dz*($Nz-1))" | bc)
for ((j=0; j<nSensorsArray; j++)); do
    for ((i=0; i<nSensorsArray; i++)); do
        xPos=$(echo "scale=4;($i*$dx*($Nx-1))/($nSensorsArray-1)" | bc)
        yPos=$(echo "scale=4;($j*$dy*($Ny-1))/($nSensorsArray-1)" | bc)
        echo "$xPos $yPos 0 $nRaysPhi $nRaysTheta -3.14 3.14 0.04 1.57" >> $INPUT_FOLDER$SENSORS
    done
done

# XZ - YZ
for ((k=1; k<nSensorsArray-1; k++)); do
    zPos=$(echo "scale=4;($k*$dz*($Nz-1))/($nSensorsArray-1)" | bc)
    # XZ
    yPos=$(echo "scale=4; 0" | bc)
    for ((i=0; i<nSensorsArray; i++)); do
        xPos=$(echo "scale=4;($i*$dx*($Nx-1))/($nSensorsArray-1)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta 0 3.14 0.04 3.1" >> $INPUT_FOLDER$SENSORS
    done
    # XY 
    xPos=$(echo "scale=4;($dx*($Nx-1))" | bc)
    for ((i=1; i<nSensorsArray-1; i++)); do
        yPos=$(echo "scale=4;($i*$dy*($Ny-1))/($nSensorsArray-1)" | bc)
        echo "0     $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.04 3.1" >> $INPUT_FOLDER$SENSORS
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta  1.57 4.71 0.04 3.1" >> $INPUT_FOLDER$SENSORS
    done 
    # XY 
    yPos=$(echo "scale=4; ($dy*($Ny-1))" | bc)
    for ((i=0; i<nSensorsArray; i++)); do
        xPos=$(echo "scale=4;($i*$dx*($Nx-1))/($nSensorsArray-1)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta 3.14 6.28 0.04 3.1" >> $INPUT_FOLDER$SENSORS
    done
done 

# XY Top
zPos=$(echo "scale=4;($dz*($Nz-1))" | bc)
for ((j=0; j<nSensorsArray; j++)); do
    for ((i=0; i<nSensorsArray; i++)); do
        xPos=$(echo "scale=4;($i*$dx*($Nx-1))/($nSensorsArray-1)" | bc)
        yPos=$(echo "scale=4;($j*$dy*($Ny-1))/($nSensorsArray-1)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta -3.14 3.14 1.57 3.1" >> $INPUT_FOLDER$SENSORS
    done
done

# Call RT solver
export OMP_NUM_THREADS=26
#RTsolver_CPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE $INPUT_FOLDER$SENSORS $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL
#RTsolver_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE $INPUT_FOLDER$SENSORS $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL

#====================
# CUDA MEMCHECK
#====================
#/opt/cuda/cuda-8.0/bin/cuda-memcheck \
RTsolver_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                                                  $INPUT_FOLDER$INITIAL_PRESSURE $INPUT_FOLDER$SENSORS \
                                                  $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL
