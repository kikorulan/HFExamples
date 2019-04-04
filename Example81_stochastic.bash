#!/bin/bash

#====================
# QSUB
#====================

#$ -P gpu
#$ -l gpu=1
#$ -l h_rt=20:00:00
#$ -l tmem=3G
#$ -N spdhg81_tau1e17
#$ -wd /home/frullan/HighFreqCode/Examples/Ex81_3D_veins_subsampled_het
#$ -S /bin/bash

# -o RTiter.txt
#$ -j y

#================================================================================
# EXAMPLE for STOCHASTIC PDHG
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
export EXAMPLE="Ex81_3D_veins_subsampled_het/"

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
export SENSORS="sensors_subsampled_3600.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_noisy5_3600sensors.dat"
export PIXEL_PRESSURE="pixelPressure_0.dat"
#export PIXEL_PRESSURE="pixelPressure_adjoint_3600sensors.dat"

# Choose GPU
export GPU_INDEX=0
# Choose mode
export MODE='-p'

#================================================================================
#=======   GRADIENT DESCENT
#================================================================================
if [ "$MODE" = "-G" ]; then
    echo "=================== GRADIENT DESCENT ===================="
    # Regularization parameters
    TAU=1e18
    LAMBDA=1e-4
    NITER=30
    # Output
    export STDOUT="stdout_GD_tau"$TAU"_lambda"$LAMBDA$"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   STOCHASTIC GRADIENT DESCENT
#================================================================================
elif [ "$MODE" = "-g" ]; then
    echo "=================== STOCHASTIC GRADIENT DESCENT ===================="
    # Regularization parameters
    TAU=4e18
    LAMBDA=1e-4
    BATCH_SIZE=100
    N_EPOCHS=30
    # Output
    export STDOUT="stdout_S-GD_tau"$TAU"_lambda"$LAMBDA"_batch"$BATCH_SIZE"_epochs"$N_EPOCHS".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $BATCH_SIZE $N_EPOCHS > $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   FISTA
#================================================================================
elif [ "$MODE" = "-F" ]; then
    echo "=================== FISTA ===================="
    # Regularization parameters
    TAU=1e18
    LAMBDA=1e-4
    NITER=30
    # Output
    export STDOUT="stdout_FISTA_tau"$TAU"_lambda"$LAMBDA"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   STOCHASTIC FISTA
#================================================================================
elif [ "$MODE" = "-f" ]; then
    echo "=================== STOCHASTIC FISTA ===================="
    # Regularization parameters
    TAU=1e18  
    LAMBDA=1e-3
    BATCH_SIZE=90
    N_EPOCHS=5
    # Output
    export STDOUT="stdout_S-FISTA_tau"$TAU"_lambda"$LAMBDA"_batch"$BATCH_SIZE"_epochs"$N_EPOCHS".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $BATCH_SIZE $N_EPOCHS> $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   PRIMAL DUAL HYBRID GRADIENT
#================================================================================
elif [ "$MODE" = "-P" ]; then
    echo "=================== PDHG ===================="
    # Regularization parameters
    SIGMA=1
    TAU=1e18
    THETA=1      
    LAMBDA=1e-4
    NITER=30
    # Output
    export STDOUT="stdout_PDHG_sigma"$SIGMA"_tau"$TAU"_theta"$THETA"_lambda"$LAMBDA"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED  \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   STOCHASTIC PRIMAL DUAL HYBRID GRADIENT
#================================================================================
elif [ "$MODE" = "-p" ]; then
    echo "=================== S-PDHG ===================="
    # Regularization parameters
    SIGMA=1
    TAU=1e17
    THETA=1    
    LAMBDA=1e-4
    BATCH_SIZE=100
    N_EPOCHS=30
    # Output
    export STDOUT="stdout_S-PDHG_sigma"$SIGMA"_tau"$TAU"_theta"$THETA"_lambda"$LAMBDA"_batch"$BATCH_SIZE"_epochs"$N_EPOCHS".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $BATCH_SIZE $N_EPOCHS > $OUTPUT_FOLDER$STDOUT
elif [ "$MODE" = "-r" ]; then
    echo "============  SINGLE FORWARD ADJOINT  ============"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE
else
    echo "Non supported mode"
fi

