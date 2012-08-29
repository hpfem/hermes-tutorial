Tutorial Part I (Linear Problems)
=================================

The first tutorial chapter begins with describing several mesh
data formats that Hermes can read, and showing how to load, refine, 
and visualize meshes. Then we learn how to set up a finite element space and 
solve a first simple problem - a Poisson equation with zero Dirichlet 
boundary conditions. After that we show how to prescribe more
general boundary conditions. We also explain how Hermes
handles numerical quadrature since this is both very important 
for higher-order finite element methods and very different 
from standard low-order FEM codes. At the end of this 
chapter we will learn how to solve systems of equations and 
axisymmetric 3D problems. 

.. toctree::
   :maxdepth: 2

   A-linear/01-mesh   
   A-linear/02-space
   A-linear/03-poisson
   A-linear/essential_and_natural_bc
   A-linear/04-bc-dirichlet
   A-linear/05-bc-neumann
   A-linear/06-bc-newton
   A-linear/quadrature
   A-linear/07-general
   A-linear/08-system
   A-linear/filters
   A-linear/09-axisym








