#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;

// This example is a matrix-free version of example 01-intro. It employs the 
// Trilinos NOX package (Trilinos needs to be installed and enabled in CMake.vars). 
//
// PDE: -div[eps_r(x,y) grad phi] = 0
//      eps_r = EPS_1 in Omega_1 (surrounding air)
//      eps_r = EPS_2 in Omega_2 (moving part of the motor)
//
// BC: phi = 0 V on Gamma_1 (left edge and also the rest of the outer boundary
//     phi = VOLTAGE on Gamma_2 (boundary of stator)
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;             // Set to "true" to enable VTK output.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.2;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H;        // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                  // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK; // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                                      // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
                                                  
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
  CustomWeakFormPoisson wf("Motor", EPS_MOTOR, "Air", EPS_AIR);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential_out("Outer", 0.0);
  DefaultEssentialBCConst<double> bc_essential_stator("Stator", VOLTAGE);
  EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_essential_out, &bc_essential_stator));

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);

  // Initialize coarse and reference mesh solution.
  Solution<double> sln, ref_sln;

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

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

    // Construct globally refined reference mesh and setup reference space.
    Space<double>* ref_space = Space<double>::construct_refined_space(&space);
    int ndof_ref = ref_space->get_num_dofs();

    // Initialize reference problem.
    info("Solving on reference mesh.");
    DiscreteProblem<double> dp(&wf, ref_space);
    
    NewtonSolverNOX<double> newton_nox(&dp);
    newton_nox.set_verbose_output(false);

    // Allocate initial coefficient vector for the Newton's method
    // on the (new) fine mesh.
    double* coeff_vec = new double[ndof_ref];

    // In the first step it is set to zero, in the following steps
    // one takes the OG projection of the previous fine mesh solution.
    if (as == 1) memset(coeff_vec, 0, ndof_ref * sizeof(double));
    else 
    {
      OGProjectionNOX<double>::project_global(ref_space, &ref_sln, coeff_vec);
    }

    // Perform Newton's iteration.
    if (!newton_nox.solve(coeff_vec)) 
      error("Newton's iteration failed.");
    else
      // Translate the resulting coefficient vector into the instance of Solution.
      Solution<double>::vector_to_solution(newton_nox.get_sln_vector(), ref_space, &ref_sln);
    
    // Project the fine mesh solution on the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjectionNOX<double>::project_global(&space, &ref_sln, &sln);

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
    Adapt<double> adaptivity(&space);
    bool solutions_for_adapt = true;
    // In the following function, the Boolean parameter "solutions_for_adapt" determines whether
    // the calculated errors are intended for use with adaptivity (this may not be the case, for example,
    // when error wrt. an exact solution is calculated). The default value is solutions_for_adapt = true,
    // The last parameter "error_flags" determine whether the total and element errors are treated as
    // absolute or relative. Its default value is error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL.
    // In subsequent examples and benchmarks, these two parameters will be often used with
    // their default values, and thus they will not be present in the code explicitly.
    double err_est_rel = adaptivity.calc_err_est(&sln, &ref_sln, solutions_for_adapt) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      space.get_num_dofs(), ref_space->get_num_dofs(), err_est_rel);

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
      done = adaptivity.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of performed adaptivity steps.
      if (done == false)  
        as++;
    }
    if (space.get_num_dofs() >= NDOF_STOP) 
      done = true;

    // Clean up.
    delete [] coeff_vec;

    // FIXME: Memory leak (no meshes, spaces or solutions are deleted 
    // during adaptivity.
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(&ref_sln);

  // Wait for all views to be closed.
  Views::View::wait();

  return 0;
}

