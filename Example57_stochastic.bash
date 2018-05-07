#!/bin/bash

#================================================================================
# EXAMPLE for STOCHASTIC PDHG
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================

# Output folder
export EXAMPLE_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex57_RT_3D_5x5x5/"
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="soundSpeed_10p.dat"
export INITIAL_PRESSURE="initialPressure_3balls.dat"
export SENSORS="sensors.dat" 
export FORWARD_SIGNAL="forwardSignal.dat"
export PIXEL_PRESSURE="pixelPressure.dat"

# Generate dimensions file
Nx=128 dx=0.001
Ny=128 dy=0.001
Nz=128 dz=0.001
cat > $INPUT_FOLDER$DIMENSIONS <<EOF
$Nx $Ny $Nz
$dx $dy $dz
EOF

# Parameters
nSensorsArray=5
nRaysPhi=1000
nRaysTheta=300
dt=1e-7
tMax=1.5e-4
# Generate sensor file
cat > $INPUT_FOLDER$SENSORS<<EOF
EOF
# XY Bottom
zPos=$(echo "scale=4;($dz*($Nz-1))" | bc)
for ((j=0; j<nSensorsArray; j++)); do
    for ((i=0; i<nSensorsArray; i++)); do
        xPos=$(echo "scale=4;($i*$dx*($Nx-1))/($nSensorsArray-1)" | bc)
        yPos=$(echo "scale=4;($j*$dy*($Ny-1))/($nSensorsArray-1)" | bc)
        echo "$xPos $yPos 0 $nRaysPhi $nRaysTheta -3.14 3.14 0.04 1.57 $dt $tMax" >> $INPUT_FOLDER$SENSORS
    done
done

# XZ - YZ
for ((k=1; k<nSensorsArray-1; k++)); do
    zPos=$(echo "scale=4;($k*$dz*($Nz-1))/($nSensorsArray-1)" | bc)
    # XZ
    yPos=$(echo "scale=4; 0" | bc)
    for ((i=0; i<nSensorsArray; i++)); do
        xPos=$(echo "scale=4;($i*$dx*($Nx-1))/($nSensorsArray-1)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta 0 3.14 0.04 3.1 $dt $tMax" >> $INPUT_FOLDER$SENSORS
    done
    # XY 
    xPos=$(echo "scale=4;($dx*($Nx-1))" | bc)
    for ((i=1; i<nSensorsArray-1; i++)); do
        yPos=$(echo "scale=4;($i*$dy*($Ny-1))/($nSensorsArray-1)" | bc)
        echo "0     $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.04 3.1 $dt $tMax" >> $INPUT_FOLDER$SENSORS
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta  1.57 4.71 0.04 3.1 $dt $tMax" >> $INPUT_FOLDER$SENSORS
    done 
    # XY 
    yPos=$(echo "scale=4; ($dy*($Ny-1))" | bc)
    for ((i=0; i<nSensorsArray; i++)); do
        xPos=$(echo "scale=4;($i*$dx*($Nx-1))/($nSensorsArray-1)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta 3.14 6.28 0.04 3.1 $dt $tMax" >> $INPUT_FOLDER$SENSORS
    done
done 

# XY Top
zPos=$(echo "scale=4;($dz*($Nz-1))" | bc)
for ((j=0; j<nSensorsArray; j++)); do
    for ((i=0; i<nSensorsArray; i++)); do
        xPos=$(echo "scale=4;($i*$dx*($Nx-1))/($nSensorsArray-1)" | bc)
        yPos=$(echo "scale=4;($j*$dy*($Ny-1))/($nSensorsArray-1)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta -3.14 3.14 1.57 3.1 $dt $tMax" >> $INPUT_FOLDER$SENSORS
    done
done

# Regularization parameters
SIGMA=1e-1
TAU=1e7
THETA=1
EPOCHS=20
# Call RT solver
export OMP_NUM_THREADS=26
RTiterative_GPU $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $EPOCHS
