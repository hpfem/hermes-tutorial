Kelly-Based *h*-Adaptivity (02-kelly)
-------------------------------------

This is yet another version of the example "01-intro" that employs 
a conventional Kelly-type error estimate to perform *h*-adaptivity. 
The Kelly error estimate is based on jumps in solution gradients 
across element edges, and it does not use the solution-pair approach
that Hermes uses by default. Most notable differences compared to 
"01-intro" are:

Calculating Kelly error estimate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Instead of a reference solution, one calculates an error estimate based
on gradient jumps::

    // Calculate element errors and total error estimate.
    info("Calculating error estimate.");
    bool ignore_visited_segments = true;
    KellyTypeAdapt<double> adaptivity(&space, ignore_visited_segments, 
      USE_EPS_IN_INTERFACE_ESTIMATOR 
      ? 
      new CustomInterfaceEstimatorScalingFunction("Motor", EPS_MOTOR, "Air", EPS_AIR)
      :
      new CustomInterfaceEstimatorScalingFunction);
    
    adaptivity.add_error_estimator_surf(new 
      Hermes::Hermes2D::BasicKellyAdapt<double>::ErrorEstimatorFormKelly());
    
    if (USE_RESIDUAL_ESTIMATOR) 
    {
      adaptivity.add_error_estimator_vol(new ResidualErrorFormMotor("Motor", EPS_MOTOR));
      adaptivity.add_error_estimator_vol(new ResidualErrorFormAir("Air", EPS_AIR));
    }

Normalizing error estimate by energy norm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
This is optional::

    if (USE_EPS_IN_INTERFACE_ESTIMATOR)
      // Use normalization by energy norm.
      adaptivity.set_error_form(new EnergyErrorForm(&wf));

The rest of the adaptivity loop is as usual.
