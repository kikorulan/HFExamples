#!/bin/bash

#===============================================================================
# This script summarizes the results stored in the given folder
#===============================================================================
SOURCE_FOLDER="/scratch0/NOT_BACKED_UP/frullan/Examples/Ex85_3D_veins_subsampled/results"
#SOURCE_FOLDER="/scratch0/NOT_BACKED_UP/frullan/Examples/Ex86_3D_veins_het/results"
#SOURCE_FOLDER="/scratch0/NOT_BACKED_UP/frullan/ExperimentalData/RD07_finger3_doubleRes_randomSubsample/results"
#SOURCE_FOLDER="/scratch0/NOT_BACKED_UP/frullan/ExperimentalData/RD10_finger2_doubleRes_subsampled/results"

echo $SOURCE_FOLDER
function containsElement () {
  local e match="$1" iter=0
  shift
  for e; 
  do 
    iter=$(($iter+1))
    if [ "$e" == "$match" ]; then
        return $iter; 
    fi
  done
  return 0
}


# Function to extract iter number
function extractIter () {
    stringIter="$1"
    if [ "${stringIter:0:1}" == "s" ]; then
        Niter=$(echo "${stringIter:8}" | cut -d'.' -f 1)
        return $Niter
    fi
    Niter=$(echo "${stringIter:4}" | cut -d'.' -f 1)
    return $Niter
}
# Change folder
cd $SOURCE_FOLDER
TAB="    "
for ELEM_i in *
do
    
    # Check if element is a folder
    if test -d "$ELEM_i" ; then
        echo $ELEM_i
        # Chech subfolders
        for ELEM_j in "$ELEM_i"/*
        do
            if test -d "$ELEM_j" ; then
                echo "$TAB$(echo $ELEM_j | cut -d'/' -f 2)"
                # Declare array
                unset ARRAY_EXPERIMENTS ARRAY_ITER
                declare -a ARRAY_EXPERIMENTS ARRAY_ITER
                # List elements in subfolder
                for ELEM_k in "$ELEM_j"/*
                do
                    # Extract element name
                    NAME_FILE=$(echo $ELEM_k | cut -d'/' -f 3)
                    ROOT_NAME=$(echo $NAME_FILE | awk '{n = split($0,SPLIT_NAME,"_"); for (i = 3; i<n; i++){ printf "%s_", SPLIT_NAME[i];} print " "}')
                    # Number of iterations
                    NITER_STRING=$(echo $NAME_FILE | awk '{n = split($0,SPLIT_NAME,"_"); print SPLIT_NAME[n]}')    
                    extractIter $NITER_STRING
                    NITER=$?
                    containsElement $ROOT_NAME "${ARRAY_EXPERIMENTS[@]}"
                    INDEX=$?
                    if [ $INDEX == 0 ]; then
                        #echo "$TAB$TAB$ROOT_NAME $NITER"
                        ARRAY_EXPERIMENTS+=($ROOT_NAME)
                        ARRAY_ITER+=($NITER)
                    else
                        if [ $NITER -ge ${ARRAY_ITER[$(($INDEX-1))]} ]; then
                            ARRAY_ITER[$(($INDEX-1))]=$NITER
                        fi
                    fi
                done
                # Show elements and iterations
                INDEX=0
                for ELEM in "${ARRAY_ITER[@]}"
                do
                    echo "$TAB$TAB${ARRAY_EXPERIMENTS[$INDEX]} $ELEM"
                    (( INDEX = INDEX+1 ))
                done
            fi
        done
    fi
done
