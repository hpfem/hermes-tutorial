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

Solving 2D Problems
-------------------

.. toctree::
    :maxdepth: 1

    src/hermes2d/P01-linear
    src/hermes2d/P02-nonlinear
    src/hermes2d/P03-transient
    src/hermes2d/P04-adaptivity
    src/hermes2d/P06-fvm-and-dg
    src/hermes2d/P07-trilinos
    src/hermes2d/P08-miscellaneous

Solving 1D Problems
-------------------

Selected 1D problems are part of the repository "hermes-examples". To clone it, type::

    git clone git://github.com/hpfem/hermes-examples.git

Solving Eigenproblems
---------------------

Eigenproblems are supported in the legacy code but not in Version 1.0,
because of their Python dependency (PySparse). Implementation of a C++
eigensolver class is scheduled for the next release.
