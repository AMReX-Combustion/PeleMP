PeleMP: Multiphysics solver for the Pele code suite
---------------------------------------------------

|AMReX Badge|
|Exascale Computing Project Badge|

`PeleMP` is the multiphysics code extension for `PeleC`, `PeleLM`, and `PeleLMeX`. Currently provides models for sprays and soot.

Getting Started
~~~~~~~~~~~~~~~

Detailed instructions are provided in the `Documentation <https://amrex-combustion.github.io/PeleMP/>`_

#. Determine which `Pele` code will work with your problem. Follow the instructions listed in that repo.

#. Set the enviroment variable, `PELEMP_HOME`, and clone a copy of `PeleMP` there ::

     export PELEMP_HOME=<path_to_PeleMP>
     git clone git@github.com:AMReX-Combustion/PeleMP.git ${PELEMP_HOME}

   or ::

     export PELEMP_HOME=<path_to_PeleMP>
     git clone https://github.com/AMReX-Combustion/PeleMP.git ${PELEMP_HOME}

Dependencies
~~~~~~~~~~~~

`PeleMP` requires the `AMReX` and `PelePhysics` libraries.

.. |AMReX Badge| image:: https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22
    :target: https://amrex-codes.github.io/amrex/
.. |Exascale Computing Project Badge| image:: https://img.shields.io/badge/supported%20by-ECP-blue
    :target: https://www.exascaleproject.org/research-project/combustion-pele/
