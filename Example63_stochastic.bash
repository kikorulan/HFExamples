#!/bin/bash

#================================================================================
# EXAMPLE for STOCHASTIC PDHG
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================

# Output folder
export EXAMPLE_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex63_3D_veins/"
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER

# Mode
export MODE='-f'

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export INITIAL_PRESSURE="initial_pressure_veins.dat"
export SENSORS="sensors.dat" 
export FORWARD_SIGNAL="forwardSignal_98sensors.dat"
export PIXEL_PRESSURE="pixelPressure.dat"

# Regularization parameters - SPDHG
#SIGMA=1e13   # 1e-2
#TAU=1e-2      # 1e10
#THETA=1      # 1
#LAMBDA=1e-3  # 1e-4
#EPOCHS=20

# Regularization parameters - PDHG
#SIGMA=0.5      # 1e-2
#TAU=1e12     # 1e10
#THETA=1      # 1
#LAMBDA=5e-3  # 1e-4
#EPOCHS=20

# Regularization parameters - FISTA
LAMBDA=1e-2
LIPSCHITZ=1e-11
NITER=50

# Call RT solver
export OMP_NUM_THREADS=26
# SPDHG
#RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
#                $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $EPOCHS
# FISTA
RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $LAMBDA $LIPSCHITZ $NITER
