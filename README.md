# PeleC-MP

PeleC-MP is the multi-physics code extension for PeleC and PeleLM.

## Getting Started

To get started, first ensure the AMReX, PelePhysics, PeleC, and PeleLM codes have been cloned and set-up properly. Then, set the enviroment variable, PELEC_MP_HOMe, and clone a copy of PeleC-MP there: ::
    export PELEC_MP_HOME=<location for PeleC-MP>
    git clone git@github.com:AMReX-Combustion/pelec-mp.git ${PELEC_MP_HOME}


## Dependencies

PeleC-MP requires the AMReX and PelePhysics. Also, the PeleC and/or PeleLM libraries are necessary, depending on which cases are being run.
