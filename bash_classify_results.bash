#!/bin/bash

#===============================================================================
# This script classifies the ".dat" files in the given folders according to their names
#===============================================================================
declare -a SOURCE_FOLDER
# List of folders
SOURCE_FOLDER[0]="/scratch0/NOT_BACKED_UP/frullan/Examples/Ex85_3D_veins_subsampled/results"
SOURCE_FOLDER[1]="/scratch0/NOT_BACKED_UP/frullan/Examples/Ex86_3D_veins_het/results"
SOURCE_FOLDER[2]="/scratch0/NOT_BACKED_UP/frullan/ExperimentalData/RD07_finger3_doubleRes_randomSubsample/results"
SOURCE_FOLDER[3]="/scratch0/NOT_BACKED_UP/frullan/ExperimentalData/RD10_finger2_doubleRes_subsampled/results"

for FOLDER in ${SOURCE_FOLDER[@]}
do
    echo $FOLDER
    # Change folder
    cd $FOLDER
    # Loop over files
    for ELEM in *
    do
        # Check if element is a file
        if test -f "$ELEM"
        then
            # Check if file has ".dat" extension
            if [ "${ELEM:(-4)}" == ".dat" ]; then
                # Extract direction (forward/adjoint) and algorithm
                DIRE=$(echo $ELEM | cut -d'_' -f 1)            
                ALGO=$(echo $ELEM | cut -d'_' -f 2)
                FLAG_MOVE=1
                # Assign folder to move
                case $DIRE in
                    forwardSignal)
                        DIRE_FOLDER="forward"
                        ;;
                    pixelPressure)
                        DIRE_FOLDER="adjoint"
                        ;;
                    *)
                        FLAG_MOVE=0
                        ;;
                esac
                # Assign subfolder to move
                case $ALGO in
                    GD)
                        ALGO_FOLDER="FB"
                        ;;
                    S-GD)
                        ALGO_FOLDER="SFB"
                        ;;
                    FISTA)
                        ALGO_FOLDER="AFB"
                        ;;
                    PDHG)
                        ALGO_FOLDER="PDHG"
                        ;;
                    S-PDHG)
                        ALGO_FOLDER="SPDHG"
                        ;;
                    *)
                        FLAG_MOVE=0
                esac
                # Move file     
                if [ $FLAG_MOVE == 1 ]; then 
                    mv $ELEM $FOLDER"/"$DIRE_FOLDER"/"$ALGO_FOLDER"/"
                    #echo $ELEM $FOLDER"/"$DIRE_FOLDER"/"$ALGO_FOLDER"/"
                fi
            fi
        fi
    done
done
