#!/bin/bash

# Tells script to stop running as soon as an error occurs
set -e
#EXEC="./PeleC2d.llvm.TPROF.MPI.ex"
EXEC="./PeleC2d.gnu.TPROF.MPI.ex"
TPD="output_files"
# For name of grid input file
# Determines box sizes on refined levels
gridsizes=(32 128)
# Box sizes for coarse level
bsize=(32 128)
mkdir -p ${TPD}

GRIDLOCS="two_d_gridfiles"
# Number of iterations, should have 4 digits
NUM_ITER_1=0005
NUM_ITER=0010
INPUT_FILE=inputs_2d

for gsi in "${!gridsizes[@]}"; do
    gs=${gridsizes[$gsi]}
    outloc=${TPD}/${gs}_BOX
    gridfile=${GRIDLOCS}/gridfile_${gs}
    regridfile=${GRIDLOCS}/gridfile_${gs}_2
    mkdir -p $outloc
    mpiexec -np 4 ${EXEC} $INPUT_FILE \
            amr.plot_file = $outloc/plt \
            amr.check_file = $outloc/chk \
            amr.checkpoint_files_output = 1\
            amr.check_int = 10000 \
            amr.plot_int = 10000 \
            amr.blocking_factor = ${bsize[$gsi]} \
            amr.max_grid_size = ${bsize[$gsi]} \
            amr.regrid_int = -1 -1 -1 \
            max_step = ${NUM_ITER_1} \
            amr.initial_grid_file = $gridfile \
            amr.regrid_file = $regridfile
    mpiexec -np 4 ${EXEC} $INPUT_FILE \
            amr.plot_file = $outloc/plt \
            amr.restart = $outloc/chk${NUM_ITER_1} \
            amr.checkpoint_files_output = 0 \
            amr.plot_int = 10000 \
            amr.blocking_factor = ${bsize[$gsi]} \
            amr.max_grid_size = ${bsize[$gsi]} \
            amr.regrid_int = -1 -1 -1 \
            max_step = ${NUM_ITER} \
            amr.initial_grid_file = $gridfile \
            amr.regrid_file = $regridfile
done
cd ${TPD}
FCOMP_EXEC=$AMREX_HOME/Tools/Plotfile/fcompare.gnu.ex
set +e
${FCOMP_EXEC} -a -r 8.E-12 ${gridsizes[0]}_BOX/plt${NUM_ITER} ${gridsizes[1]}_BOX/plt${NUM_ITER}>comp.log
cat comp.log
