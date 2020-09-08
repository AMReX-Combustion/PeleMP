#!/bin/bash

# Tells script to stop running as soon as an error occurs
set -e
EXEC="./PeleC2d.gnu.TPROF.MPI.ex"
TPD="output_files"
# For name of grid input file
# Determines box sizes on refined levels
gridsizes=(32 64 128)
# Box sizes for coarse level
bsize=(32 64 256)
mkdir -p ${TPD}

GRIDLOCS="two_d_gridfiles"
# Number of iterations, should have 6 digits
NUM_ITER=000005
INPUT_FILE=inputs_2d

for gsi in "${!gridsizes[@]}"; do
    gs=${gridsizes[$gsi]}
    outloc=${TPD}/${gs}_BOX
    gridfile=${GRIDLOCS}/gridfile_${gs}
    mkdir -p $outloc
    mpiexec -np 4 ${EXEC} $INPUT_FILE \
            amr.plot_file = $outloc/plt \
            amr.blocking_factor = ${bsize[$gsi]} \
            amr.max_grid_size = ${bsize[$gsi]} \
            max_step = ${NUM_ITER} \
            amr.plot_int = ${NUM_ITER} \
            amr.initial_grid_file = $gridfile
done
cd ${TPD}
FCOMP_EXEC=/Users/ldowen/Codes/amrex/Tools/Plotfile/fcompare.gnu.ex
${FCOMP_EXEC} -a ${gridsizes[0]}_BOX/plt${NUM_ITER} ${gridsizes[1]}_BOX/plt${NUM_ITER}>out1.log
${FCOMP_EXEC} -a ${gridsizes[0]}_BOX/plt${NUM_ITER} ${gridsizes[2]}_BOX/plt${NUM_ITER}>out2.log
${FCOMP_EXEC} -a ${gridsizes[2]}_BOX/plt${NUM_ITER} ${gridsizes[1]}_BOX/plt${NUM_ITER}>out3.log
cat out1.log
cat out2.log
cat out3.log
