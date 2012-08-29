Trilinos - Linear (01-trilinos-linear)
--------------------------------------

The purpose of this example is to show how to use Trilinos for linear PDE problems where
NOX only does one nonlinear iteration. Thus this is a good example to test various preconditioners
for the linear problem. First we use the Newton's method in Hermes (assembling via the DiscreteProblem 
class and matrix problem solution via UMFpack). Second, assembling is done using the DiscreteProblem 
class in Hermes and the discrete problem is solved using the Trilinos NOX solver (using Newton's 
method or JFNK, with or without preconditioning).


Model problem
~~~~~~~~~~~~~

The PDE solved is 

.. math::
    -\Delta u - f = 0

with an exact solution 

.. math::
    u(x,y) = x^2 + y^2.

The first part (UMFpack) needs not be discussed, let's proceed directly to NOX: 

Calculate initial vector for NOX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As NOX is iterative, the initial condition matters::

    // Calculate initial vector for NOX.
    info("Projecting to obtain initial vector for the Newton's method.");
    ZeroSolution init_sln(&mesh);
    OGProjection<double> ogProjection; ogProjection.project_global(&space, &init_sln, coeff_vec);

Initializing NOX
~~~~~~~~~~~~~~~~

::

    // Initialize the NOX solver.
    info("Initializing NOX.");
    NewtonSolverNOX<double> nox_solver(&dp2);
    nox_solver.set_output_flags(message_type);

Setting additional parameters for NOX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    // Set various NOX parameters.
    nox_solver.set_ls_type(iterative_method);
    nox_solver.set_ls_tolerance(ls_tolerance);
    nox_solver.set_conv_iters(max_iters);
    if (flag_absresid)
      nox_solver.set_conv_abs_resid(abs_resid);
    if (flag_relresid)
      nox_solver.set_conv_rel_resid(rel_resid);

Setting a preconditioner
~~~~~~~~~~~~~~~~~~~~~~~~

::

    // Choose preconditioning.
    MlPrecond<double> pc("sa");
    if (PRECOND)
    {
      if (TRILINOS_JFNK) nox_solver.set_precond(pc);
      else nox_solver.set_precond(preconditioner);
    }

See NOX documentation for more preconditioning choices.

Calling NOX
~~~~~~~~~~~

Now we are ready to call the NOX solver to assemble the discrete problem and solve it::

    // Assemble and solve using NOX.
    try
    {
      nox_solver.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("NOX failed.");
    }

Collecting output info
~~~~~~~~~~~~~~~~~~~~~~

::

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      nox_solver.get_num_iters(), nox_solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      nox_solver.get_num_lin_iters(), nox_solver.get_achieved_tol());

Sample results
~~~~~~~~~~~~~~

You should see the following result:

.. figure:: 01-trilinos-linear/1.png
   :align: center
   :scale: 75% 
   :figclass: align-center
   :alt: Sample result
