#!/bin/bash

# Output folder
export EXAMPLE_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex52_RT_3D/"
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
0     0     0     1000 300     0 1.57 0.04 1.57 1e-7 1e-4
0.127 0     0     1000 300  1.57 3.14 0.04 1.57 1e-7 1e-4
0     0.127 0     1000 300 -1.57    0 0.04 1.57 1e-7 1e-4
0.127 0.127 0     1000 300  3.14 4.71 0.04 1.57 1e-7 1e-4
0     0     0.127 1000 300     0 1.57 1.57 3.10 1e-7 1e-4
0.127 0     0.127 1000 300  1.57 3.14 1.57 3.10 1e-7 1e-4
0     0.127 0.127 1000 300 -1.57    0 1.57 3.10 1e-7 1e-4
0.127 0.127 0.127 1000 300  3.14 4.71 1.57 3.10 1e-7 1e-4
EOF

# Call RT solver
RTsolver_CPU $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE $INPUT_FOLDER$SENSORS
