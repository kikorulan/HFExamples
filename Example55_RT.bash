#!/bin/bash

#================================================================================
# EXAMPLE 53
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================

# Output folder
export EXAMPLE_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex55_RT_3D_half/"
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"03/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="soundSpeed_half.dat"
export INITIAL_PRESSURE="initialPressure_2balls.dat"
export SENSORS="sensors.dat" 
export FORWARD_SIGNAL="forwardSignal.dat"

# Mode
export MODE="-a"

# Generate dimensions file
Nx=128 dx=0.001
Ny=128 dy=0.001
Nz=128 dz=0.001
cat > $INPUT_FOLDER$DIMENSIONS <<EOF
$Nx $Ny $Nz
$dx $dy $dz
EOF

# Parameters
nSensorsArray=3
nRaysPhi_aux=1000
nRaysTheta_aux=300
nRaysPhi=1000
nRaysTheta=300
dt=1e-7
tMax=1.5e-4
# Generate sensor file
cat > $INPUT_FOLDER$SENSORS<<EOF
EOF

# XZ
for ((k=0; k<nSensorsArray; k++)); do
    zPos=$(echo "scale=4;($k*$dz*($Nz-1))/($nSensorsArray-1)" | bc)
    # XZ
    yPos=$(echo "scale=4; 0" | bc)
    for ((i=0; i<nSensorsArray; i++)); do
        xPos=$(echo "scale=4;($i*$dx*($Nx-1))/($nSensorsArray-1)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta 0 3.14 0.04 3.1 $dt $tMax" >> $INPUT_FOLDER$SENSORS
    done
    # XZ
    yPos=$(echo "scale=4; ($dy*($Ny-1))" | bc)
    for ((i=0; i<nSensorsArray; i++)); do
        xPos=$(echo "scale=4;($i*$dx*($Nx-1))/($nSensorsArray-1)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta 3.14 6.28 0.04 3.1 $dt $tMax" >> $INPUT_FOLDER$SENSORS
    done
done 


# Call RT solver
export OMP_NUM_THREADS=26
RTsolver_CPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE $INPUT_FOLDER$SENSORS $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL
