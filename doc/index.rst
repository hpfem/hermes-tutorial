===============
Hermes Tutorial
===============

.. raw:: latex
 
   \floatplacement{figure}{H}

Thank you for your interest in Hermes!

Hermes is a C++ library for rapid development of adaptive *hp*-FEM and *hp*-DG solvers,
with emphasis on nonlinear, time-dependent, multi-physics problems.

This document is organized as follows: 

* Section 1 provides general information about Hermes and the computational methods it uses,
  and gives references to underlying scientific articles.
* Section 2 describes how to install Hermes on various hardware platforms, and how to 
  install matrix solver packages and various optional packages. 
* Section 3 explains how to use Git and Github, and how you can contribute to the project if interested.
* Section 4 contains a tutorial to Hermes2D. Please read this tutorial first even if you are 
  interested in 1D or 3D problems, since the syntax is virtually the same. The tutorial 
  will walk you in small steps through the solution
  of linear, nonlinear, and time-dependent problems from various engineering and scientific areas, 
  using higher-order elements and adaptivity algorithms, and solving multiphysics coupled problems. 
* Section 5 shows how to solve 1D problems.

This document is under continuous development. If you find bugs, typos, dead links 
and such, please report them to the 
`Hermes2D mailing list <http://groups.google.com/group/hermes2d/>`_ -- thanks!

Introduction
------------

.. toctree::
    :maxdepth: 1

    src/about-hermes
    src/math-background
    src/web-access
    src/citing-hermes

Installation
------------

.. toctree::
    :maxdepth: 1

    src/installation/linux
    src/installation/mac
    src/installation/win-cygwin
    src/installation/win-msvc
    src/installation/matrix_solvers
    src/installation/cython_installation
    src/installation/exodusII_netcdf

Collaboration
-------------

.. toctree::
    :maxdepth: 1

    src/collaboration
    src/editing_sphinx

Tutorial Examples
-----------------

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

There is no separate Hermes1D library -- one-dimensional problems are solved using Hermes2D. 

.. toctree::
    :maxdepth: 1

    src/hermes1d/examples.rst
    src/hermes1d/quantum-notes.rst

