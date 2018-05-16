#!/bin/bash

#================================================================================
# EXAMPLE for STOCHASTIC PDHG
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================

# Output folder
export EXAMPLE_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex60_3D_veins/"
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export INITIAL_PRESSURE="initial_pressure_4balls.dat"
export SENSORS="sensors.dat" 
export FORWARD_SIGNAL="forwardSignal_56sensors.dat"
export PIXEL_PRESSURE="pixelPressure.dat"

# Regularization parameters
SIGMA=1e-1
TAU=1e11
THETA=1
LAMBDA=1e-2
EPOCHS=20
# Call RT solver
export OMP_NUM_THREADS=26
RTiterative_GPU $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $EPOCHS
