#!/bin/bash

#================================================================================
# EXAMPLE for STOCHASTIC PDHG
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================

export EXAMPLE="Ex72_3D_veins_heterogeneous/"

# Output folder
if [ "$HOSTNAME" = "maryam.cs.ucl.ac.uk" ]; then
    export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/"
elif [ "$HOSTNAME" = "hannover" ]; then
    export HOST_FOLDER="/home/wontek/sharedWK/Examples/"
elif [ "$HOSTNAME" = "miller.local" ] || [ "$HOSTNAME" = "armstrong.local" ]; then
    export HOST_FOLDER="/home/frullan/HighFreqCode/Examples/"
fi
export EXAMPLE_FOLDER=$HOST_FOLDER$EXAMPLE
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
export PIXEL_PRESSURE="pressure_adjoint_kWave_57600sensors.dat"
export STDOUT="stdout"$1".txt"

# Choose GPU
export GPU_INDEX=0

#================================================================================
#=======   GRADIENT DESCENT
#================================================================================
if [ "$1" = "-G" ]; then
    echo "=================== GRADIENT DESCENT ===================="
    # Regularization parameters
    TAU=5e18
    LAMBDA=1e-2
    NITER=50
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT

#================================================================================
#=======   STOCHASTIC GRADIENT DESCENT
#================================================================================
elif [ "$1" = "-g" ]; then
    echo "=================== STOCHASTIC GRADIENT DESCENT ===================="
    # Regularization parameters
    TAU=5e18
    LAMBDA=3e-4
    NITER=50
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT

#================================================================================
#=======   FISTA
#================================================================================
elif [ "$1" = "-F" ]; then
    echo "=================== FISTA ===================="
    # Regularization parameters
    TAU=5e18
    LAMBDA=1e-2
    NITER=50
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT

#================================================================================
#=======   STOCHASTIC FISTA
#================================================================================
elif [ "$1" = "-f" ]; then
    echo "=================== STOCHASTIC FISTA ===================="
    # Regularization parameters
    TAU=5e18   
    LAMBDA=3e-3
    NITER=50
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT

#================================================================================
#=======   PRIMAL DUAL HYBRID GRADIENT
#================================================================================
elif [ "$1" = "-P" ]; then
    echo "=================== PDHG ===================="
    # Regularization parameters
    SIGMA=1e0
    TAU=5e18       # 1e10
    THETA=1        # 1
    LAMBDA=1e-2    # 1e-4
    NITER=50
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   STOCHASTIC PRIMAL DUAL HYBRID GRADIENT
#================================================================================
elif [ "$1" = "-p" ]; then
    echo "=================== S-PDHG ===================="
    # Regularization parameters
    SIGMA=1e0
    TAU=1e18
    THETA=1    
    LAMBDA=1e-4
    EPOCHS=20
    # RUN
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $EPOCHS > $OUTPUT_FOLDER$STDOUT

elif [ "$1" = "-r" ]; then
    echo "============  SINGLE FORWARD ADJOINT  ============"
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE
else
    echo "Non supported mode"
fi

    
