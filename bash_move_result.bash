#!/bin/bash

#===============================================================================
# This script moves a specified result
#===============================================================================
ROOT_FOLDER="/scratch0/NOT_BACKED_UP/frullan/"
# Source
SOURCE_FOLDER="Examples/"
SOURCE_SUBFOLDER="Ex85_3D_veins_subsampled/"
#SOURCE_SUBFOLDER="Ex86_3D_veins_het/"
SOURCE_RESULT_FOLDER="results/"
# Destination
DEST_FOLDER="Examples/"
#DEST_SUBFOLDER="Ex85_3D_veins_subsampled/"
DEST_SUBFOLDER="Ex86_3D_veins_het/"
#DEST_FOLDER="ExperimentalData/"
#DEST_SUBFOLDER="RD07_finger3_doubleRes_randomSubsample/"
#DEST_SUBFOLDER="RD10_finger2_doubleRes_subsampled/"
DEST_RESULT_FOLDER="results/"


# Parameters
LAMBDA='1e-4'
TAU='1.6e2'
SIGMA='1'
THETA='1'
BATCH='100'

# Algorithm
ALGO='GD'

#========== REMOVE
# Change folder
cd $ROOT_FOLDER

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
mv "$SOURCE_FOLDER""$SOURCE_SUBFOLDER""$SOURCE_RESULT_FOLDER"forward/"$ALGO_FOLDER"/forwardSignal_"$NAME_FILE"* "$DEST_FOLDER""$DEST_SUBFOLDER""$DEST_RESULT_FOLDER"forward/"$ALGO_FOLDER"/
mv "$SOURCE_FOLDER""$SOURCE_SUBFOLDER""$SOURCE_RESULT_FOLDER"adjoint/"$ALGO_FOLDER"/pixelPressure_"$NAME_FILE"* "$DEST_FOLDER""$DEST_SUBFOLDER""$DEST_RESULT_FOLDER"adjoint/"$ALGO_FOLDER"/
#rm -f adjoint/"$ALGO_FOLDER"/pixelPressure_"$NAME_FILE"*
