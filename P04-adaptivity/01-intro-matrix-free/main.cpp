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
                                                  
// Problem parameters.
const double EPS0 = 8.863e-12;
const double VOLTAGE = 50.0;
const double EPS_MOTOR = 10.0 * EPS0;
const double EPS_AIR = 1.0 * EPS0;

// NOX parameters.
const bool TRILINOS_JFNK = true;                  // true = Jacobian-free method (for NOX),
                                                  // false = Newton (for NOX).
const bool PRECOND = true;                        // Preconditioning by jacobian in case of JFNK (for NOX),
                                                  // default ML preconditioner in case of Newton.
const char* iterative_method = "GMRES";           // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "AztecOO";           // Name of the preconditioner employed by AztecOO 
                                                  // Possibilities: None" - No preconditioning. 
                                                  // "AztecOO" - AztecOO internal preconditioner.
                                                  // "New Ifpack" - Ifpack internal preconditioner.
                                                  // "ML" - Multi level preconditione
unsigned message_type = NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration | NOX::Utils::InnerIteration | NOX::Utils::Parameters | NOX::Utils::LinearSolverDetails;
                                                  // NOX error messages, see NOX_Utils.h.

double ls_tolerance = 1e-5;                       // Tolerance for linear system.
unsigned flag_absresid = 0;                       // Flag for absolute value of the residuum.
double abs_resid = 1.0e-3;                        // Tolerance for absolute value of the residuum.
unsigned flag_relresid = 1;                       // Flag for relative value of the residuum.
double rel_resid = 1.0e-2;                        // Tolerance for relative value of the residuum.
int max_iters = 100;                              // Max number of iterations.


int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Motor", EPS_MOTOR, "Air", EPS_AIR, TRILINOS_JFNK);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential_out("Outer", 0.0);
  DefaultEssentialBCConst<double> bc_essential_stator("Stator", VOLTAGE);
  EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_essential_out, &bc_essential_stator));

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);

  // Initialize coarse and fine mesh solution.
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
  Space<double>* ref_space_new = NULL, *ref_space_prev = NULL;
  do
  {
    info("---- Adaptivity step %d:", as);
    
    // Time measurement.
    cpu_time.tick();

    // Backup fine mesh space if it already exists.
    if (ref_space_new != NULL) 
    {
      if (ref_space_prev != NULL) delete ref_space_prev;
      ref_space_prev = ref_space_new;
    }

    // Construct (new) fine mesh and setup (new) fine mesh space.
    ref_space_new = Space<double>::construct_refined_space(&space);
    int ndof_ref = ref_space_new->get_num_dofs();

    // Initialize (new) fine mesh problem.
    DiscreteProblem<double> dp(&wf, ref_space_new);
    
    // Allocate initial coefficient vector for the Newton's method
    // on the (new) fine mesh.
    double* coeff_vec = new double[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(double));

    // Initialize the NOX solver with the vector "coeff_vec".
    info("Initializing NOX.");
    NewtonSolverNOX<double> newton_nox(&dp);
    newton_nox.set_verbose_output(true);
    newton_nox.set_output_flags(message_type);
    newton_nox.set_ls_type(iterative_method);
    newton_nox.set_ls_tolerance(ls_tolerance);
    newton_nox.set_conv_iters(max_iters);
    if (flag_absresid)
      newton_nox.set_conv_abs_resid(abs_resid);
    if (flag_relresid)
      newton_nox.set_conv_rel_resid(rel_resid);

    // Transfer previous fine mesh solution to new fine mesh.
    if (as > 1)
    {
      info("Transferring previous fine mesh solution to new fine mesh.");
      OGProjectionNOX<double>::project_global(ref_space_new, &ref_sln, coeff_vec);
    }

    // Choose preconditioning.
    MlPrecond<double> pc("sa");
    if (PRECOND)
    {
      if (TRILINOS_JFNK) newton_nox.set_precond(pc);
      else newton_nox.set_precond(preconditioner);
    }

    // Perform Newton's iteration.
    info("Performing JFNK on new fine mesh.");
    if (!newton_nox.solve(coeff_vec)) 
      error("JFNK iteration failed.");
    else
    {
      // Translate the resulting coefficient vector into the instance of Solution.
      Solution<double>::vector_to_solution(newton_nox.get_sln_vector(), ref_space_new, &ref_sln);

      // Output.
      info("Number of nonlin iterations: %d (norm of residual: %g)", 
        newton_nox.get_num_iters(), newton_nox.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
        newton_nox.get_num_lin_iters(), newton_nox.get_achieved_tol());
    }

    // Project the fine mesh solution on the coarse mesh.
    info("Projecting fine mesh solution to coarse mesh.");
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
      space.get_num_dofs(), ref_space_new->get_num_dofs(), err_est_rel);

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

  // Show the last fine mesh solution - final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(&ref_sln);

  // Wait for all views to be closed.
  Views::View::wait();

  return 0;
}

