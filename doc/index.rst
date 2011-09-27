===============
Hermes Tutorial
===============

.. raw:: latex
 
   \floatplacement{figure}{H}

Thank you for your interest in Hermes!

Hermes is a C++ library for rapid development of adaptive *hp*-FEM and *hp*-DG solvers,
with emphasis on nonlinear, time-dependent, multi-physics problems.

This document contains a beginner's tutorial. It will walk you in small steps through the solution
of linear, nonlinear, and time-dependent problems from various engineering and scientific areas, 
using higher-order elements and adaptivity algorithms, and solving multiphysics coupled problems. 
Section 2 shows how to solve 1D problems.

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
    src/hermes2d/P05-eigenproblems
    src/hermes2d/P06-fvm-and-dg
    src/hermes2d/P07-trilinos
    src/hermes2d/P08-miscellaneous

Solving 1D Problems
-------------------

Selected 1D problems are part of the repository "hermes-examples". To clone it, type::

    git clone http://git.hpfem.org/git/hermes-examples.git

The repository contains a folder doc/ with Sphinx documentation. To build the docs,
type "make html" in that directory. 
