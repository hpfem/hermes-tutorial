Adaptive Multimesh *hp*-FEM Example (03-system)
-----------------------------------------------

Model problem
~~~~~~~~~~~~~

We consider a simplified version of the Fitzhugh-Nagumo equation.
This equation is a prominent example of activator-inhibitor systems in two-component reaction-diffusion 
equations, It describes a prototype of an excitable system (e.g., a neuron) and its stationary form 
is

.. math::

    -d^2_u \Delta u - f(u) + \sigma v - g_1 = 0,\\
    -d^2_v \Delta v - u + v - g_2 = 0.

Here the unknowns $u, v$ are the voltage and $v$-gate, respectively.
The nonlinear function 

.. math::

    f(u) = \lambda u - u^3 - \kappa
 
describes how an action potential travels through a nerve. Obviously this system is nonlinear.

Exact solution
~~~~~~~~~~~~~~

In order to make it simpler for this tutorial, we replace the function $f(u)$ with just $u$:

.. math::

    f(u) = u.

Our computational domain is the square $(-1,1)^2$ and we consider zero Dirichlet conditions 
for both $u$ and $v$. In order to enable fair convergence comparisons, we will use the following 
functions as the exact solution:

.. math::

    u(x,y) = \cos\left(\frac{\pi}{2}x\right) \cos\left(\frac{\pi}{2}y\right),\\
    v(x,y) = \hat u(x) \hat u(y)

where

.. math::

    \hat u(x) = 1 - \frac{e^{kx} + e^{-kx}}{e^k + e^{-k}}

is the exact solution of the one-dimensional singularly perturbed 
problem 

.. math::

    -u'' + k^2 u - k^2 = 0

in $(-1,1)$, equipped with zero Dirichlet boundary conditions. 

The following two figures show the solutions $u$ and $v$. Notice their 
large qualitative differences: While $u$ is smooth in the entire domain, 
$v$ has a thin boundary layer along the boundary:

.. figure:: 03-system/solution_u.png
   :align: center
   :scale: 40% 
   :figclass: align-center
   :alt: Solution

.. figure:: 03-system/solution_v.png
   :align: center
   :scale: 40% 
   :figclass: align-center
   :alt: Solution

Manufactured right-hand side
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The source functions $g_1$ and $g_2$ are obtained by inserting $u$ and $v$ 
into the PDE system. These functions are not extremely pretty, but they 
are not too bad either - see files definitions.h and definitions.cpp.

Weak forms
~~~~~~~~~~

The weak forms can be found in the files definitions.h and definitions.cpp.
Beware that although each of the forms is actually symmetric, one cannot use the 
HERMES_SYM flag as in the elasticity equations, since it has a different 
meaning.

Adaptivity loop
~~~~~~~~~~~~~~~

The adaptivity workflow is standard. First we construct the reference spaces::

    // Construct globally refined reference mesh and setup reference space.
    Hermes::vector<Space<double> *>* ref_spaces = 
      Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&u_space, &v_space));

Next we initialize the discrete problem on the fine meshes::

    // Initialize reference problem.
    DiscreteProblem<double> dp(&wf, *ref_spaces);

Then we initialize the Newton solver::

    NewtonSolver<double> newton(&dp, matrix_solver);
    newton.set_verbose_output(false);

And we solve the problem using the Newton's method::

    // Perform Newton's iteration.
    try
    {
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    }

Next we translate the coefficient vector into the two Solutions,
project reference solutions to the coarse meshes, calculate error 
estimates, calculate exact errors (optional), and 
adapt the coarse meshes. For details on the last steps see
the filemain.cpp.


Sample results
~~~~~~~~~~~~~~

Now we can show some numerical results. 
First let us show the resulting meshes for $u$ and $v$ obtained using 
conventional (single-mesh) hp-FEM: **9,330 DOF** (4665 for each solution component). 

.. figure:: 03-system/mesh_single.png
   :align: center
   :scale: 40% 
   :figclass: align-center
   :alt: Mesh

.. raw:: html

   <hr style="clear: both; visibility: hidden;">

.. figure:: 03-system/mesh_single.png
   :align: center
   :scale: 40% 
   :figclass: align-center
   :alt: Mesh

.. raw:: html

   <hr style="clear: both; visibility: hidden;">

Next we show the resulting meshes for $u$ and $v$ obtained using 
the multimesh hp-FEM: **1,723 DOF** (49 DOF for $u$ and $1,673$ for $v$). 

.. figure:: 03-system/mesh_multi_u.png
   :align: center
   :scale: 40% 
   :figclass: align-center
   :alt: Mesh

.. raw:: html

   <hr style="clear: both; visibility: hidden;">

.. figure:: 03-system/mesh_multi_v.png
   :align: center
   :scale: 40% 
   :figclass: align-center
   :alt: Mesh

.. raw:: html

   <hr style="clear: both; visibility: hidden;">

Finally let us compare the DOF and CPU convergence graphs 
for both cases:

DOF convergence graphs:

.. figure:: 03-system/conv_dof.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. figure:: 03-system/conv_cpu.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: CPU convergence graph.
