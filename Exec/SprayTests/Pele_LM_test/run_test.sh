#!/bin/sh

# Tells script to stop running as soon as an error occurs
set -e
#EXEC="./PeleLM2d.gnu.MPI.ex"
EXEC="./PeleLM2d.llvm.MPI.ex"
TPD="output_files"
# For name of grid input file
# Determines box sizes on refined levels
gridsizes=(32 128)
# Box sizes for coarse level
bsize=(32 128)
mkdir -p ${TPD}

GRIDLOCS="two_d_gridfiles"
# Number of iterations, should have 6 digits
NUM_ITER_1=00005
NUM_ITER=00010
#INPUT_FILE=closed-input
INPUT_FILE=open-input

for gsi in "${!gridsizes[@]}"; do
    gs=${gridsizes[$gsi]}
    outloc=${TPD}/${gs}_BOX
    gridfile=${GRIDLOCS}/gridfile_${gs}
    regridfile=${gridfile}_2
    mkdir -p $outloc
    # Run initial
    mpiexec -np 4 ${EXEC} $INPUT_FILE \
            amr.plot_file = $outloc/plt \
            amr.check_file = $outloc/chk \
            amr.checkpoint_files_output = 1 \
            amr.blocking_factor = ${bsize[$gsi]} \
            amr.max_grid_size = ${bsize[$gsi]} \
            max_step = ${NUM_ITER_1} \
            amr.initial_grid_file = $gridfile \
            amr.regrid_file = $regridfile \
            amr.regrid_int = 4 4 4 4
    # Run from restart
    mpiexec -np 4 ${EXEC} $INPUT_FILE \
            amr.plot_file = $outloc/plt \
            amr.check_file = $outloc/chk \
            amr.restart = $outloc/chk${NUM_ITER_1} \
            amr.blocking_factor = ${bsize[$gsi]} \
            amr.max_grid_size = ${bsize[$gsi]} \
            max_step = ${NUM_ITER} \
            amr.initial_grid_file = $gridfile \
            amr.regrid_file = $regridfile \
            amr.regrid_int = 4 4 4 4
done
cd ${TPD}
# Must fill in the path to fcompare here
FCOMP_EXEC=$AMREX_HOME/Tools/Plotfile/fcompare.gnu.ex
set +e
${FCOMP_EXEC} -a -r 8.E-12 ${gridsizes[0]}_BOX/plt${NUM_ITER} ${gridsizes[1]}_BOX/plt${NUM_ITER}>comp.log
cat comp.log
