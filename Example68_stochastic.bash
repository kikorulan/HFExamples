#!/bin/bash

#================================================================================
# EXAMPLE for STOCHASTIC PDHG
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================

# Output folder
export EXAMPLE_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex68_3D_veins_resize/"
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER

# Mode
export MODE='-p'

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export INITIAL_PRESSURE="initial_pressure_veins_80x240x240.dat"
export SENSORS="sensors_subsampled.dat" 
export FORWARD_SIGNAL="forwardSignal_reference.dat"
export PIXEL_PRESSURE="pixelPressure.dat"

if [ "$MODE" = "-s" ]; then
    echo "=================== SPDHG ===================="
    # Regularization parameters - SPDHG
    SIGMA=1e13   # 1e-2
    TAU=1e-2      # 1e10
    THETA=1      # 1
    LAMBDA=1e-3  # 1e-4
    EPOCHS=20
    # RUN
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $EPOCHS
elif [ "$MODE" = "-p" ]; then
    echo "=================== PDHG ===================="
    # Regularization parameters - PDHG
    SIGMA=0.5      # 1e-2
    TAU=1e14       # 1e10
    THETA=1        # 1
    LAMBDA=5e-3    # 1e-4
    NITER=20
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $NITER
elif [ "$MODE" = "-f" ]; then
    echo "=================== FISTA ===================="
    # Regularization parameters - FISTA
    LAMBDA=1e-2
    LIPSCHITZ=1e-11
    NITER=50
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $LAMBDA $LIPSCHITZ $NITER
fi
