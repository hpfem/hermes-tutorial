Using NURBS in Computation (02-nurbs)
-------------------------------------

Model problem
~~~~~~~~~~~~~

Solved is a simple Poisson equation with constant right-hand
side and homogeneous Dirichlet boundary conditions.

The domain is a rectangle (0,2) x (0, 1) whose upper
edge is a NURBS curve. There are three mesh files
in this example: domain-1.mesh (one control point),
domain-2.mesh (two control points), and domain-3.mesh
(three control points). One of these files needs to be 
selected on line 15 in main.cpp::

    const char* mesh_file = "domain-1.mesh";          // One control point.
    //const char* mesh_file = "domain-2.mesh";          // Two control points.
    //const char* mesh_file = "domain-3.mesh";          // Three control points.

Example 1 (One inner control point)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Snippet from the file domain-1.mesh::

    degree = 2                # Degree should be equal to the 
    num_inner_points = 1      # number of inner points plus one.
                              
    inner_points =
    {
      { 1.5, 2.0, 2.0 }       # x, y, weight
    } 
    knots = 
    {
      0, 0, 0, 1, 1, 1        
    }
    curves =
    {
      {2, 3, degree, inner_points, knots} 
    }

Result:

.. figure:: 02-nurbs/1.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: NURBS with one control point.

Example 2 (Two inner control points)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Snippet from the file domain-3.mesh::

    degree = 3
    num_inner_points = 2
    inner_points =
    {
      { 1.5, 1.5, 1.0 },
      { 0.5, 0.5, 1.0 }
    } 
    knots = 
    {
      0, 0, 0, 1, 1, 1
    }
    curves =
    {
      {2, 3, degree, inner_points, knots} 
    }

Result:

.. figure:: 02-nurbs/2.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: NURBS with two control points.


Example 3 (Three inner control points)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Snippet from the file domain-2.mesh::

    degree = 4
    num_inner_points = 3
    inner_points =
    {
      { 1.5, 1.5, 1.0 },
      { 1.0, -1.0, 1.0 },
      { 0.5, 1.5, 1.0 }
    } 
    knots = 
    {
      0, 0, 0, 1, 1, 1
    }
    curves =
    {
      {2, 3, degree, inner_points, knots} 
    }

Result:

.. figure:: 02-nurbs/3.png
   :align: center
   :scale: 45% 
   :figclass: align-center
   :alt: NURBS with three control points.




