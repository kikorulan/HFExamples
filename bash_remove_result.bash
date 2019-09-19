#!/bin/bash

#===============================================================================
# This script removes the specified result
#===============================================================================
ROOT_FOLDER="/scratch0/NOT_BACKED_UP/frullan/"
EXAMPLE_FOLDER="Examples/"
EXAMPLE_SUBFOLDER="Ex85_3D_veins_subsampled/"
#EXAMPLE_SUBFOLDER="Ex86_3D_veins_het/"
#EXAMPLE_FOLDER="ExperimentalData/"
#EXAMPLE_SUBFOLDER="RD10_finger2_doubleRes_subsampled/"
RESULT_FOLDER="results/"

# Parameters
LAMBDA='1e-3'
TAU='6.4e2'
SIGMA='2'
THETA='1'
BATCH='100'

# Algorithm
ALGO='S-GD'

#========== REMOVE
# Change folder
cd $ROOT_FOLDER$EXAMPLE_FOLDER$EXAMPLE_SUBFOLDER$RESULT_FOLDER

# Switch between algorithms
case $ALGO in
    GD) 
        ALGO_FOLDER="FB"
        NAME_FILE="GD_tau"$TAU"_lambda"$LAMBDA"_iter"
        ;;
    S-GD)
        ALGO_FOLDER="SFB"
        NAME_FILE="S-GD_tau"$TAU"_lambda"$LAMBDA"_batch"$BATCH"_subepoch"
        ;;
    FISTA)
        ALGO_FOLDER="AFB"
        NAME_FILE="FISTA_tau"$TAU"_lambda"$LAMBDA"_iter"
        ;;
    PDHG)
        ALGO_FOLDER="PDHG"
        NAME_FILE="PDHG_sigma"$SIGMA"_tau"$TAU"_theta"$THETA"_lambda"$LAMBDA"_iter"
        ;;
    S-PDHG)
        ALGO_FOLDER="SPDHG"
        NAME_FILE="S-PDHG_sigma"$SIGMA"_tau"$TAU"_theta"$THETA"_lambda"$LAMBDA"_batch"$BATCH"_subepoch"
        ;;
    *)
        echo "Unknown algorithm"
        exit 0
esac


# Remove from forward and adjoint
rm -f forward/"$ALGO_FOLDER"/forwardSignal_"$NAME_FILE"*
rm -f adjoint/"$ALGO_FOLDER"/pixelPressure_"$NAME_FILE"*
