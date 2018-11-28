#!/bin/bash

#====================
# QSUB
#====================

#$ -l tmem=3G
#$ -l h_vmem=3G
#$ -l h_rt=3:00:00
#$ -l gpu=1

# OPTIONAL

#$ -S /bin/bash
#$ -j y
#$ -N RTsolver

#================================================================================
# EXAMPLE 68
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
#The code you want to run now goes here.
export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
export EXAMPLE="Ex72_3D_veins_heterogeneous/"
# Output folder
if [ "$HOSTNAME" = "maryam.cs.ucl.ac.uk" ]; then
    export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/"
elif [ "$HOSTNAME" = "hannover" ]; then
    export HOST_FOLDER="/home/wontek/sharedWK/Examples/"
else
    export HOST_FOLDER="/home/frullan/HighFreqCode/Examples/"
fi
export EXAMPLE_FOLDER=$HOST_FOLDER$EXAMPLE
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER

# Choose GPU
export GPU_INDEX=0

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export INITIAL_PRESSURE="initial_pressure_veins_80x240x240.dat"
export SENSORS="sensors_subsampled_1600.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_1600sensors.dat"
export STDOUT="stdout-adjoint-1600sensors.txt"

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

#==============================
# SENSORS
#==============================
nSensorsArray=40
nRaysPhi=1024 
nRaysTheta=1024
dt=1.3e-8
tMax=8.0836e-06
# Generate sensor file
cat > $INPUT_FOLDER$SENSORS<<EOF
EOF
# Step and tMax
echo "$dt $tMax 0 0 0 0 0 0 0" >> $INPUT_FOLDER$SENSORS
# YZ
for ((k=0; k<nSensorsArray; k++)); do
    zPos=$(echo "scale=4;($k*$dz*($Nz-1))/($nSensorsArray-1)" | bc)
    for ((i=0; i<nSensorsArray; i++)); do
        yPos=$(echo "scale=4;($i*$dy*($Ny-1))/($nSensorsArray-1)" | bc)
        echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.04 3.1" >> $INPUT_FOLDER$SENSORS
    done 
done 

#====================
# RUN 
#====================
RTsolver_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
             $INPUT_FOLDER$SENSORS $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL  > $OUTPUT_FOLDER$STDOUT
