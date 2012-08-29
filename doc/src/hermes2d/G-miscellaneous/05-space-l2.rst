Space L2 (05-space-l2)
----------------------

The L2 space is used, for example, in `Navier-Stokes equations where 
it keeps the velocity discreetely divergence-free. This example shows how to 
create an L2 space, visualize finite element basis functions, and perform 
an orthogonal L2-projection of a continuous function onto the FE space.

First, an L2 space is created as the reader expects::

    // Create an L2 space with default shapeset.
    L2Space<double> space(&mesh, P_INIT);

The function to be projected is::

    CustomExactSolution sln_exact(&mesh);

See formula in the file definitions.cpp. The projection is done as follows::

    OGProjection<double>::project_global(&space, &sln_exact, &sln, matrix_solver);

Sample basis functions visualized using the BaseView class:

.. figure:: 05-space-l2/fn0.png
   :align: center
   :scale: 45% 
   :figclass: align-center
   :alt: Sample basis function

.. figure:: 05-space-l2/fn1.png
   :align: center
   :scale: 45% 
   :figclass: align-center
   :alt: Sample basis function

.. figure:: 05-space-l2/fn2.png
   :align: center
   :scale: 45% 
   :figclass: align-center
   :alt: Sample basis function

.. figure:: 05-space-l2/fn3.png
   :align: center
   :scale: 45% 
   :figclass: align-center
   :alt: Sample basis function

The projection (note that this is a discontinuous function):

.. figure:: 05-space-l2/sol.png
   :align: center
   :scale: 45% 
   :figclass: align-center
   :alt: Projection
