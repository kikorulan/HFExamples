#!/bin/bash

#====================
# QSUB
#====================

#$ -P gpu
#$ -l gpu=1
#$ -l h_rt=10:00:00
#$ -l tmem=3G
#$ -N spdhg81_tau2e19_sigma1
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

export EXAMPLE="Ex82_TVanalysis/"

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

# Choose GPU
export GPU_INDEX=1

# Parameters
export N_ITER_TV=20
export N_ITER_OUT=3

#====================
# LAMBDA = 1E-3
#====================
export LAMBDA=1e-4
PIXEL_PRESSURE_IN="u0.dat"
for ((i=0; i<N_ITER_OUT; i++)); do
    iNext=$(echo "$i+1" | bc)
    echo $iNext
    PIXEL_PRESSURE_OUT="u"$iNext"_lambda"$LAMBDA".dat"
    TVapp_GPU $INPUT_FOLDER$DIMENSIONS $OUTPUT_FOLDER$PIXEL_PRESSURE_IN $LAMBDA $N_ITER_TV $OUTPUT_FOLDER$PIXEL_PRESSURE_OUT
    PIXEL_PRESSURE_IN=$PIXEL_PRESSURE_OUT
done
