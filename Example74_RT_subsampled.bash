#!/bin/bash

#====================
# QSUB
#====================

#$ -P gpu
#$ -l gpu=1
#$ -l h_rt=10:50:0
#$ -l tmem=3G
#$ -N RTsolver
#$ -wd /home/frullan/C++
#$ -S /bin/bash

#$ -o RTsolver.txt
#$ -j y

#================================================================================
# EXAMPLE 74
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
#The code you want to run now goes here.
export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
export EXAMPLE="Ex74_3D_thinveins/"
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


# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export INITIAL_PRESSURE="initial_pressure_veins_80x240x240.dat"
export SENSORS="sensors_subsampled_14400.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_14400sensors.dat"
export STDOUT="stdout-forward.txt"

# Mode
export MODE="-f"
export GPU_INDEX=0
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
nSensorsArray=120
nRaysPhi=1024 
nRaysTheta=1024
dt=1.6667e-8
tMax=8.0836e-06
# Generate sensor file
cat > $INPUT_FOLDER$SENSORS<<EOF
EOF
# Step and tMax
echo "$dt $tMax 0 0 0 0 0 0 0" >> $INPUT_FOLDER$SENSORS
# YZ
for ((k=0; k<nSensorsArray; k++)); do
    zPos=$(echo "scale=6;($k*$dz*($Nz-1))/($nSensorsArray-1)" | bc)
    for ((i=0; i<nSensorsArray; i++)); do
        yPos=$(echo "scale=6;($i*$dy*($Ny-1))/($nSensorsArray-1)" | bc)
        echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.04 3.1" >> $INPUT_FOLDER$SENSORS
    done 
done 

#====================
# RUN 
#====================
RTsolver_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
             $INPUT_FOLDER$SENSORS $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL > $OUTPUT_FOLDER$STDOUT
