#!/bin/bash

#===============================================================================
# This script copies the ".dat" files in a folder and 
# moves them to a specified DESTINATION folder
#===============================================================================
SOURCE_FOLDER="/home/frullan/HighFreqCode/Examples/Ex85_3D_veins_subsampled/output_data/"
DEST_FOLDER="/scratch0/NOT_BACKED_UP/frullan/Examples/Ex85_3D_veins_subsampled/results/"
SLEEP_TIME=5
FILE_LIST="ListOfFiles.txt"

# Change folder
cd $SOURCE_FOLDER
# Create log of copied files
rm -f $FILE_LIST
# Infinite loop
while true
do
    # Loop over files
    for ELEM in *
    do
        # Check if element is a file
        if test -f "$ELEM"
        then
            # Check if file has ".dat" extension
            if [ "${ELEM:(-4)}" == ".dat" ]; then
                echo $ELEM >> $FILE_LIST
                scp $ELEM frullan@ember:$DEST_FOLDER
                rm $ELEM
            fi
        fi
    done

    # Sleep
    sleep "$SLEEP_TIME"h
done