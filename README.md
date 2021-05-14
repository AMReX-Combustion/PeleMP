# PeleMP

PeleMP is the multi-physics code extension for PeleC and PeleLM.

## Getting Started

To get started, first ensure the AMReX, IAMR, PelePhysics, PeleC, and PeleLM codes have been cloned and set-up properly. Then, set the enviroment variable, PELEMP_HOME, and clone a copy of PeleMP there
```
    export PELEMP_HOME=<location for PeleMP>
    git clone git@github.com:AMReX-Combustion/pelemp.git ${PELEMP_HOME}
```

## Dependencies

PeleMP requires the AMReX and PelePhysics libraries. Also, the PeleC and/or PeleLM libraries are necessary, depending on which cases are being run.
