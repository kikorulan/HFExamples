#!/bin/bash

#===============================================================================
# This script copies the ".dat" files in a folder and 
# moves them to a specified DESTINATION folder
#===============================================================================
declare -a SOURCE_FOLDER
N_FOLDERS=4
# Source folders
SOURCE_FOLDER[0]="/home/frullan/HighFreqCode/Examples/Ex85_3D_veins_subsampled/output_data/"
SOURCE_FOLDER[1]="/home/frullan/HighFreqCode/Examples/Ex86_3D_veins_het/output_data/"
SOURCE_FOLDER[2]="/home/frullan/HighFreqCode/ExperimentalData/RD07_finger3_doubleRes_randomSubsample/output_data/"
SOURCE_FOLDER[3]="/home/frullan/HighFreqCode/ExperimentalData/RD10_finger2_doubleRes_subsampled/output_data/"

# Destination folders
DEST_FOLDER[0]="/scratch0/NOT_BACKED_UP/frullan/Examples/Ex85_3D_veins_subsampled/results/"
DEST_FOLDER[1]="/scratch0/NOT_BACKED_UP/frullan/Examples/Ex86_3D_veins_het/results/"
DEST_FOLDER[2]="/scratch0/NOT_BACKED_UP/frullan/ExperimentalData/RD07_finger3_doubleRes_randomSubsample/results/"
DEST_FOLDER[3]="/scratch0/NOT_BACKED_UP/frullan/ExperimentalData/RD10_finger2_doubleRes_subsampled/results/"
SLEEP_TIME=1
FILE_LIST="ListOfFiles.txt"

# Create log of copied files
rm -f $FILE_LIST
# Infinite loop
while true
do
    for ((i=0;i<N_FOLDERS;i++)); do
        # Change folder
        cd ${SOURCE_FOLDER[$i]}
        # Loop over files
        for ELEM in *
        do
            # Check if element is a file
            if test -f "$ELEM"
            then
                # Check if file has ".dat" extension
                if [ "${ELEM:(-4)}" == ".dat" ]; then
                    echo $ELEM >> $FILE_LIST
                    scp $ELEM frullan@ember:${DEST_FOLDER[$i]}
                    rm $ELEM
                fi
            fi
        done
    done
    # Sleep
    sleep "$SLEEP_TIME"h
done
