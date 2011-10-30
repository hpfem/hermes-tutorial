#define HERMES_REPORT_ALL
#include "definitions.h"

// This example shows how to run adaptive h-FEM driven by the Kelly estimator and
// set its basic control parameters. The underlying problem is the same as in 
// P04-adaptivity/01-intro, i.e., a planar model of an electrostatic micromotor (MEMS).
//
// In this example, adaptive mesh refinement is not based on comparing solutions on 
// coarse and refined meshes as in previous tutorial, but rather on evaluating an
// error estimate for each element. The basic error estimator developed by Kelly and
// co-workers for elliptic problems is represented by class BasicKellyAdapt (its 
// base class KellyTypeAdapt may be used to adjust the estimator to more general
// problems). It calculates the error of an element as L2 norm of jumps of 
// solution accross element edges. In order to fully conform to the theory, it may 
// also optionally include the L2 norm of the equation residual (variable 
// USE_RESIDUAL_ESTIMATOR).
//
// Elements with highest error estimates will be refined according to the selected 
// strategy (variables THRESHOLD and STRATEGY). See the User Documentation for more 
// details. 
//
// Uniform initial polynomial degree of mesh elements can be set using the variable 
// P_INIT. Additional control parameters are available, these will be demonstrated
// in the following tutorial examples. 
//
// In this example, two types of convergence graphs are created -- error estimate 
// wrt. the number of degrees of freedom (DOF), and error estimate wrt. CPU time. 
// Later we will show how to output the error wrt. exact solution when exact
// solution is available.
//
// PDE: -div[eps_r(x,y) grad phi] = 0
//      eps_r = EPS_1 in Omega_1 (surrounding air)
//      eps_r = EPS_2 in Omega_2 (moving part of the motor)
//
// BC: phi = 0 V on Gamma_1 (left edge and also the rest of the outer boundary
//     phi = VOLTAGE on Gamma_2 (boundary of stator)
//
// The following parameters can be changed:

// Set to "false" to suppress Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;           
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;             
// Initial polynomial degree of mesh elements.
const int P_INIT = 2;                             
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;                     
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
const int STRATEGY = 0;                           
// Add also the norm of residual to the error estimate of each element.
const bool USE_RESIDUAL_ESTIMATOR = false;          
// If true, the interface estimator is defined by jumps of fluxes; otherwise
// jumps of normal derivatives are used.
const bool USE_EPS_IN_INTERFACE_ESTIMATOR = true;   
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// Stopping criterion for adaptivity.
const double ERR_STOP = 5e-6;                     
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK; 
                                                  
// Problem parameters.
const double EPS0 = 8.863e-12;
const double VOLTAGE = 50.0;
const double EPS_MOTOR = 10.0 * EPS0;
const double EPS_AIR = 1.0 * EPS0;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("motor.mesh", &mesh);

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Motor", EPS_MOTOR, "Air", EPS_AIR, &mesh);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential_out("Outer", 0.0);
  DefaultEssentialBCConst<double> bc_essential_stator("Stator", VOLTAGE);
  EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_essential_out, &bc_essential_stator));

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);

  // Initialize coarse and reference mesh solution.
  Solution<double> sln;

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 410, 600));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  Views::OrderView  oview("Polynomial orders", new Views::WinGeom(420, 0, 400, 600));

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  // Time measurement.
  TimePeriod cpu_time;

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    
    // Time measurement.
    cpu_time.tick();

    // Initialize reference problem.
    info("Solving.");
    DiscreteProblem<double> dp(&wf, &space);
    
    NewtonSolver<double> newton(&dp, matrix_solver);
    newton.set_verbose_output(false);

    // Initial coefficient vector for the Newton's method.  
    int ndof = space.get_num_dofs();
    double* coeff_vec = new double[ndof];
    memset(coeff_vec, 0, ndof * sizeof(double));

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
    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);
    
    // Time measurement.
    cpu_time.tick();

    // VTK output.
    if (VTK_VISUALIZATION) 
    {
      // Output solution in VTK format.
      Views::Linearizer lin;
      char* title = new char[100];
      sprintf(title, "sln-%d.vtk", as);
      lin.save_solution_vtk(&sln, title, "Potential", false);
      info("Solution in VTK format saved to file %s.", title);

      // Output mesh and element orders in VTK format.
      Views::Orderizer ord;
      sprintf(title, "ord-%d.vtk", as);
      ord.save_orders_vtk(&space, title);
      info("Element orders in VTK format saved to file %s.", title);
    }

    // View the coarse mesh solution and polynomial orders.
    if (HERMES_VISUALIZATION) 
    {
      sview.show(&sln);
      oview.show(&space);
    }

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);
    
    // Calculate element errors and total error estimate.
    info("Calculating error estimate.");
    bool ignore_visited_segments = true;
    KellyTypeAdapt<double> adaptivity(&space, ignore_visited_segments, 
                                      USE_EPS_IN_INTERFACE_ESTIMATOR 
                                        ? 
                                          new CustomInterfaceEstimatorScalingFunction("Motor", EPS_MOTOR, "Air", EPS_AIR)
                                        :
                                          new CustomInterfaceEstimatorScalingFunction);
    
    adaptivity.add_error_estimator_surf(new Hermes::Hermes2D::BasicKellyAdapt<double>::ErrorEstimatorFormKelly());
    
    if (USE_RESIDUAL_ESTIMATOR) 
    {
      adaptivity.add_error_estimator_vol(new ResidualErrorFormMotor("Motor", EPS_MOTOR));
      adaptivity.add_error_estimator_vol(new ResidualErrorFormAir("Air", EPS_AIR));
    }
    
    if (USE_EPS_IN_INTERFACE_ESTIMATOR)
      // Use normalization by energy norm.
      adaptivity.set_error_form(new EnergyErrorForm(&wf));
    
    // Note that there is only one solution (the only one available) passed to BasicKellyAdapt::calc_err_est
    // and there is also no "solutions_for_adapt" parameter. The last parameter, "error_flags", is left 
    // intact and has the meaning explained in P04-adaptivity/01-intro. You may however still call
    // adaptivity.calc_err_est with a solution, an exact solution and solutions_for_adapt=false to calculate
    // error wrt. an exact solution (if provided). 
    double err_est_rel = adaptivity.calc_err_est(&sln, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

    // Report results.
    info("ndof: %d, err_est_rel: %g%%", space.get_num_dofs(), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    cpu_time.tick();    
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    graph_dof.add_values(space.get_num_dofs(), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    
    // Skip the time spent to save the convergence graphs.
    cpu_time.tick(HERMES_SKIP);

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) 
      done = true;
    else
    {
      info("Adapting coarse mesh.");
      
      // h-refinement is automatically selected here, no need for parameter "selector".
      done = adaptivity.adapt(THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of performed adaptivity steps.
      if (done == false)  
        as++;
    }
    if (space.get_num_dofs() >= NDOF_STOP) 
      done = true;

    // Clean up.
    delete [] coeff_vec;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  // The final result has already been shown in the final step of the adaptivity loop, so we only
  // adjust the title and hide the mesh here.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);

  // Wait for all views to be closed.
  Views::View::wait();

  return 0;
}

