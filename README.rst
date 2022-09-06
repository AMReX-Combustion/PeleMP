PeleMP
------

*The multiphysics solver for the Pele code suite*

`PeleMP` is the multiphysics code extension for `PeleC`, `PeleLM`, and `PeleLMeX`.

Getting Started
~~~~~~~~~~~~~~~

1. Determine which `Pele` code will work with your problem. Follow the instructions listed in that repo.

2. Set the enviroment variable, `PELEMP_HOME`, and clone a copy of `PeleMP` there
::
   export PELEMP_HOME=<path_to_PeleMP>
   git clone git@github.com:AMReX-Combustion/PeleMP.git ${PELEMP_HOME}

or ::

  export PELEMP_HOME=<path_to_PeleMP>
  git clone https://github.com/AMReX-Combustion/PeleMP.git ${PELEMP_HOME}

Dependencies
~~~~~~~~~~~~

`PeleMP` requires the `AMReX` and `PelePhysics` libraries.
