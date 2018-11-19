#!/bin/bash

#================================================================================
# EXAMPLE for STOCHASTIC PDHG
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================

# Output folder
if [ "$HOSTNAME" = "maryam.cs.ucl.ac.uk" ]; then
    export EXAMPLE_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex68_3D_veins_resize/"
elif [ "$HOSTNAME" = "hannover" ]; then
    export EXAMPLE_FOLDER="/home/wontek/sharedWK/Examples/Ex68_3D_veins_resize/"
fi
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER

# Mode provided in the first argument of the script

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export INITIAL_PRESSURE="initial_pressure_veins_80x240x240.dat"
export SENSORS="sensors_subsampled_1600.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_1600sensors.dat"
export PIXEL_PRESSURE="pressure_kWave_adjoint_57600sensors.dat"
#export PIXEL_PRESSURE="pixelPressure_PDHG_40iter.dat"


#================================================================================
#=======   GRADIENT DESCENT
#================================================================================
if [ "$1" = "-G" ]; then
    echo "=================== GRADIENT DESCENT ===================="
    # Regularization parameters - WORKS
    TAU=1e18
    LAMBDA=1e-2 # 1e-2
    NITER=50
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER

#================================================================================
#=======   STOCHASTIC GRADIENT DESCENT
#================================================================================
elif [ "$1" = "-g" ]; then
    echo "=================== STOCHASTIC GRADIENT DESCENT ===================="
    # Regularization parameters
    TAU=1e18
    LAMBDA=4e-2 #3.5e-2
    NITER=50
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER

#================================================================================
#=======   FISTA
#================================================================================
elif [ "$1" = "-F" ]; then
    echo "=================== FISTA ===================="
    # Regularization parameters
    TAU=1e18
    LAMBDA=1e-2
    NITER=50
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER

#================================================================================
#=======   STOCHASTIC FISTA
#================================================================================
elif [ "$1" = "-f" ]; then
    echo "=================== STOCHASTIC FISTA ===================="
    # Regularization parameters
    TAU=2e18    # 1e12
    LAMBDA=4e-2 # 1e-2
    NITER=50
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER

#================================================================================
#=======   PRIMAL DUAL HYBRID GRADIENT
#================================================================================
elif [ "$1" = "-P" ]; then
    echo "=================== PDHG ===================="
    # Regularization parameters
    SIGMA=1e-20      # 1e-2
    TAU=1e-20        # 1e10
    THETA=1        # 1
    LAMBDA=1e-17    # 1e-4
    NITER=50
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $NITER
#================================================================================
#=======   STOCHASTIC PRIMAL DUAL HYBRID GRADIENT
#================================================================================
elif [ "$1" = "-p" ]; then
    echo "=================== SPDHG ===================="
    # Regularization parameters
    SIGMA=1e13   # 1e-2
    TAU=1e-2      # 1e10
    THETA=1      # 1
    LAMBDA=1e-3  # 1e-4
    EPOCHS=20
    # RUN
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $EPOCHS

elif [ "$1" = "-r" ]; then
    echo "============  SINGLE FORWARD ADJOINT  ============"
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE
else
    echo "Non supported mode"
fi

    
