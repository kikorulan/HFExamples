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


# Assign files
export DIMENSIONS="dimensions.dat"
export CUBE="pixelPressure.dat"
# Parameters
export LAMBDA=1e-14
export NITER=3000
export OUTPUT="pixelPressure_TVdenoised_$LAMBDA-$NITER.dat"


TVapp_GPU $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$CUBE $LAMBDA $NITER $OUTPUT_FOLDER$OUTPUT

    
