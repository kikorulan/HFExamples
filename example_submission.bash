#   This is the most basic QSUB file needed for this cluster.
#   Further examples can be found under /share/apps/examples
#   Most software is NOT in your PATH but under /share/apps
#
#   For further info please read http://hpc.cs.ucl.ac.uk
#
#   NOTE hash dollar is a scheduler diredctive not a comment.


# These are flags you must include - Two memory and one runtime.
# Runtime is either seconds or hours:min:sec

#$ -l tmem=2G
#$ -l h_vmem=2G
#$ -l h_rt=3:00:00
#$ -l gpu=1

#These are optional flags but you problably want them in all jobs

#$ -S /bin/bash
#$ -j y
#$ -N RTsolver

#The code you want to run now goes here.
export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"

bash /home/frullan/
hostname
date
