Using the Implicit Euler Method (01-implicit-euler)
---------------------------------------------------

Model problem
~~~~~~~~~~~~~

This section describes how a time-integration method can be built into the weak 
formulation of a PDE (usually this is called *Rothe's method*). We will illustrate this on 
the implicit Euler method applied to a heat transfer equation. The equation describes in a naive approximation 
how the St. Vitus cathedral in Prague responds to changes in the surrounding 
air temperature during one 24-hour cycle. The geometry is shown below:

.. figure:: 01-implicit-euler/mesh.png
   :align: center
   :scale: 30% 
   :figclass: align-center
   :alt: Model geometry.

We assume the standard heat transfer equation

.. math::
    :label: eqvit1

       c \varrho\frac{\partial T}{\partial t} - \lambda \Delta T = 0

equipped with a Dirichlet condition

.. math::

     T = T_{init}

on the bottom edge $\Gamma_{ground}$ and a Newton condition

.. math::

     \frac{\partial T}{\partial \nu} = \alpha(T_{ext}(t) - T)

on the rest of the boundary $\Gamma_{air}$. Here, $c$ is the heat capacity of the material,
$\varrho$ the material density, $\lambda$ the thermal conductivity,
$T_{init}$ the fixed temperature on the
ground (same as the initial temperature of the building), and $\alpha$
the heat transfer coefficient 
between the building and the surrounding air. The surrounding air temperature
$T_{ext}$ is time-dependent of the form

.. math::

     T_{ext}(t) = T_{init} + 10\sin(2\pi t/T_{final}),

where $T_{final}$ is 24 hours (translated into seconds).

Equation :eq:`eqvit1` is equipped with an initial condition of the
form

.. math::

     T(x,y,0) = T_{init}(x,y) \ \ \ \mbox{in} \ \Omega.

Replacing in :eq:`eqvit1` the temporal derivative with the backward time difference, 
we obtain

.. math::

     c \varrho\frac{T^{n+1} - T^n}{\tau} - \lambda \Delta T^{n+1} = 0.

Weak formulation
~~~~~~~~~~~~~~~~

The corresponding weak formulation reads

.. math::

     \int_{\Omega} c \varrho\frac{T^{n+1} - T^n}{\tau}v + \int_{\Omega} \lambda \nabla T^{n+1}\cdot \nabla v + \int_{\Gamma_{air}} \alpha \lambda T^{n+1}v - \int_{\Gamma_{air}} \alpha \lambda T_{ext}(t^{n+1})v = 0.

Defining weak forms
~~~~~~~~~~~~~~~~~~~

The weak formulation is a combination of default and custom weak forms, see files
definitions.h and definitions.cpp.

Passing and accessing previous time level solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Notice how previous time level solution is passed into the volumetric vector form::

    CustomVectorFormVol* vec_form_vol = new CustomVectorFormVol(0, time_step);
    vec_form_vol->ext.push_back(prev_time_sln);
    add_vector_form(vec_form_vol);

and also how it is accessed from inside the weak form::

    Func<double> *temp_prev_time = ext->fn[0];

Keeping the Jacobian if it did not change
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As this problem is linear, the Jacobian matrix just needs to be constructed once
at the beginning, and it will not change during the computation. Thus instead of 
the usual method NewtonSolver::solve() we use the method NewtonSolver::solve_keep_jacobian()::

    // This is important in time-dependent examples, these methods are shared by all the 'calculation' classes, i.e. RungeKutta, NewtonSolver, PicardSolver, LinearSolver, DiscreteProblem, DiscreteProblemLinear.
			
    // Perform Newton's iteration.
    try
    {
      newton.solve_keep_jacobian(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      
    }

Sample results
~~~~~~~~~~~~~~

Sample temperature distribution is shown below: 

.. figure:: 01-implicit-euler/vitus1.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: sample result


