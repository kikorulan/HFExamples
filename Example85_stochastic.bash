#!/bin/bash

#====================
# QSUB
#====================

# #$ -P gpu
#$ -l gpu=1
#$ -l h_rt=200:00:00
#$ -l tmem=3G
#$ -N spdhg85_tau1.6e2_sigma1e-1
#$ -wd /home/frullan/HighFreqCode/Examples/Ex85_3D_veins_subsampled
#$ -S /bin/bash

# -o RTiter.txt
#$ -j y

#================================================================================
# EXAMPLE 85 ITERATIVE
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
#The code you want to run now goes here.
if [ "$HOSTNAME" = "miller.local" ] || [ "$HOSTNAME" = "armstrong.local" ]; then
    export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build_miller/bin:$PATH"
else
    export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
fi

export EXAMPLE="Ex85_3D_veins_subsampled/"

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
export SENSORS="sensors_subsampled_3600.dat"
export FORWARD_SIGNAL="forwardSignal_reference_noisy5_3600sensors.dat"
export PIXEL_PRESSURE="pixelPressure_0.dat"

# Machine
echo $HOSTNAME
# Choose GPU
export GPU_INDEX=0
# Choose mode
export MODE='-p'

# Parameters
SIGMA=1e-1
TAU=1.6e2
THETA=1    
LAMBDA=1e-4
BATCH_SIZE=100
NITER=10

#=======   GRADIENT DESCENT
if [ "$MODE" = "-G" ]; then
    echo "=================== GRADIENT DESCENT ===================="
    # Output
    export STDOUT="stdout_GD_tau"$TAU"_lambda"$LAMBDA$"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#=======   STOCHASTIC GRADIENT DESCENT
elif [ "$MODE" = "-g" ]; then
    echo "=================== STOCHASTIC GRADIENT DESCENT ===================="
    # Output
    export STDOUT="stdout_S-GD_tau"$TAU"_lambda"$LAMBDA"_batch"$BATCH_SIZE"_epochs"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $BATCH_SIZE $NITER > $OUTPUT_FOLDER$STDOUT
#=======   FISTA
elif [ "$MODE" = "-F" ]; then
    echo "=================== FISTA ===================="
    # Output
    export STDOUT="stdout_FISTA_tau"$TAU"_lambda"$LAMBDA"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#=======   PRIMAL DUAL HYBRID GRADIENT
elif [ "$MODE" = "-P" ]; then
    echo "=================== PDHG ===================="
    # Output
    export STDOUT="stdout_PDHG_sigma"$SIGMA"_tau"$TAU"_theta"$THETA"_lambda"$LAMBDA"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#=======   STOCHASTIC PRIMAL DUAL HYBRID GRADIENT
elif [ "$MODE" = "-p" ]; then
    echo "=================== S-PDHG ===================="
    # Output
    export STDOUT="stdout_S-PDHG_sigma"$SIGMA"_tau"$TAU"_theta"$THETA"_lambda"$LAMBDA"_batch"$BATCH_SIZE"_epochs"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $BATCH_SIZE $NITER > $OUTPUT_FOLDER$STDOUT
elif [ "$MODE" = "-r" ]; then
    echo "============  SINGLE FORWARD ADJOINT  ============"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE
else
    echo "Non supported mode"
fi

    
