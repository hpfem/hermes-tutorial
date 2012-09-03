===============
Hermes Tutorial
===============

.. raw:: latex
 
   \floatplacement{figure}{H}

Thank you for your interest in Hermes!

This document contains a beginner's tutorial. It will walk you in small steps through the solution
of linear, nonlinear, and time-dependent problems from various engineering and scientific areas, 
using higher-order elements and adaptivity algorithms, and solving multiphysics coupled problems. 

The document is under continuous development. If you find bugs, typos, dead links 
and such, please report them to the 
`Hermes2D mailing list <http://groups.google.com/group/hermes2d/>`_.
    
Installation
--------------------
- installation is very simple and has the following steps:

| 1) clone the repository
| 2) on Windows (and on other platforms as well if necessary), copy CMake.vars.example to CMake.vars and adjust it.
| 3) HERMES_DIRECTORY and HERMES_INCLUDE_PATH should point to lib and include directories where you installed HERMES.
|    typically, these will be "lib" and "include" subdirectories of the directory specified by TARGET_ROOT in your library repository.
| 4) the DEP_ROOT variable should in most cases point to the same variable in CMake.vars in your library repository.
| 5) the same for the PTHREAD_ROOT variable
| 6!!) important step is to set WITH_TRILINOS to the same value as you set it in your library repository. Otherwise, linking issues will occur.

Solving 2D Problems
-------------------

.. toctree::
    :maxdepth: 1

    src/hermes2d/A-linear
    src/hermes2d/B-nonlinear
    src/hermes2d/C-transient
    src/hermes2d/D-adaptivity
    src/hermes2d/E-fvm-and-dg
    src/hermes2d/F-trilinos
    src/hermes2d/G-miscellaneous

Solving 1D Problems
-------------------

Selected 1D problems are part of the repository "hermes-examples". To clone it, type::

    git clone git://github.com/hpfem/hermes-examples.git

Solving Eigenproblems
---------------------

Eigenproblems are supported in the legacy code but not in the OpenMP version,
because of their Python dependency (PySparse). Implementation of a C++
eigensolver class is scheduled for near future.
