PeleMP: Multiphysics solver for the Pele code suite
---------------------------------------------------

|AMReX Badge|
|Exascale Computing Project Badge|

`PeleMP` was the multiphysics code extension for `PeleC`, `PeleLM`, and `PeleLMeX`. It provided models for sprays and soot.

WARNING
-------

The physics modules and documentation of PeleMP have been moved to
`PelePhysics <https://github.com/AMReX-Combustion/PelePhysics>`_,
and this repository is now archived. Further development is continuing
within PelePhysics. The test cases have moved to the PeleC and PeleLMeX repositorties.

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

Citation
~~~~~~~~

To cite the Soot and Spray capabilities from PeleMP, please use the following `Journal of Fluids Engineering article <https://doi.org/10.1115/1.4064494>`_: ::

  @article{owen2023pelemp,
    title={PeleMP: The Multiphysics Solver for the Combustion Pele Adaptive Mesh Refinement Code Suite},
    author={Owen, Landon D and Ge, Wenjun and Rieth, Martin and Arienti, Marco and Esclapez, Lucas and S Soriano, Bruno and Mueller, Michael E and Day, Marc and Sankaran, Ramanan and Chen, Jacqueline H},
    journal={Journal of Fluids Engineering},
    pages={1--41},
    year={2023}
  }


Dependencies
~~~~~~~~~~~~

`PeleMP` requires the `AMReX` and `PelePhysics` libraries.

.. |AMReX Badge| image:: https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22
    :target: https://amrex-codes.github.io/amrex/
.. |Exascale Computing Project Badge| image:: https://img.shields.io/badge/supported%20by-ECP-blue
    :target: https://www.exascaleproject.org/research-project/combustion-pele/
