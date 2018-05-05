#!/bin/bash

# Output folder
export EXAMPLE_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex42_RT_3D/"
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="soundSpeed_10p.dat"
export INITIAL_PRESSURE="initialPressure_3balls.dat"
export SENSORS="sensors.dat" 

# Generate dimensions file
cat > $INPUT_FOLDER$DIMENSIONS <<EOF
128 128 128
0.001 0.001 0.001
EOF

# Generate sensors file
cat > $INPUT_FOLDER$SENSORS <<EOF
0 0.063 0.063 1000 300 -3.14 3.14 0.1 3.1 1e-7 1e-4
EOF

# Call RT solver
RTsolver_CPU $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE $INPUT_FOLDER$SENSORS
