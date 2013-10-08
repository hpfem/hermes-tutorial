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

// Set to "false" to suppress Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;           
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;             
// Initial polynomial degree of mesh elements.
const int P_INIT = 2;                             
// Parameter influencing the candidate selection.
const double THRESHOLD = 0.2;                      
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
// H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO_H;                 
// Stopping criterion for adaptivity.
const double ERR_STOP = 1.0;   

// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Selector.
H1ProjBasedSelector<double> selector(CAND_LIST);

// Problem parameters.
const double EPS0 = 8.863e-12;
const double VOLTAGE = 50.0;
const double EPS_MOTOR = 10.0 * EPS0;
const double EPS_AIR = 1.0 * EPS0;

// NOX parameters.
// true = Jacobian-free method (for NOX),
// false = Newton (for NOX).
const bool TRILINOS_JFNK = true;                  
// Preconditioning by jacobian in case of JFNK (for NOX),
// default ML preconditioner in case of Newton.
const bool PRECOND = true;                        
// Name of the iterative method employed by AztecOO (ignored
// by the other solvers). 
// Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* iterative_method = "GMRES";           
// Name of the preconditioner employed by AztecOO 
// Possibilities: None" - No preconditioning. 
// "AztecOO" - AztecOO internal preconditioner.
// "New Ifpack" - Ifpack internal preconditioner.
// "ML" - Multi level preconditione
const char* preconditioner = "AztecOO";           
// NOX error messages, see NOX_Utils.h.
unsigned message_type = NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration | NOX::Utils::InnerIteration | NOX::Utils::Parameters | NOX::Utils::LinearSolverDetails;
                                                  
// Tolerance for linear system.
double ls_tolerance = 1e-5;                       
// Flag for absolute value of the residuum.
unsigned flag_absresid = 0;                       
// Tolerance for absolute value of the residuum.
double abs_resid = 1.0e-3;                        
// Flag for relative value of the residuum.
unsigned flag_relresid = 1;                       
// Tolerance for relative value of the residuum.
double rel_resid = 1.0e-2;                        
// Max number of iterations.
int max_iters = 100;                              

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Motor", EPS_MOTOR, "Air", EPS_AIR, TRILINOS_JFNK);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential_out("Outer", 0.0);
  DefaultEssentialBCConst<double> bc_essential_stator("Stator", VOLTAGE);
  EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_essential_out, &bc_essential_stator));

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Initialize coarse and fine mesh solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>), ref_sln(new Solution<double>);

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
  Hermes::Mixins::TimeMeasurable cpu_time;

  // Adaptivity loop:
  int as = 1; bool done = false;
  Space<double>* ref_space_new = NULL, *ref_space_prev = NULL;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);
    
    // Time measurement.
    cpu_time.tick();

    // Backup fine mesh space if it already exists.
    if (ref_space_new != NULL) 
    {
      if (ref_space_prev != NULL) delete ref_space_prev;
      ref_space_prev = ref_space_new;
    }

    // Construct (new) fine mesh and setup (new) fine mesh space.
    ref_space_new = Space<double>::construct_refined_space(space);
    int ndof_ref = ref_space_new->get_num_dofs();

    // Initialize (new) fine mesh problem.
    DiscreteProblem<double> dp(&wf, ref_space_new);
    
    // Allocate initial coefficient vector for the Newton's method
    // on the (new) fine mesh.
    double* coeff_vec = new double[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(double));

    // Initialize the NOX solver with the vector "coeff_vec".
    Hermes::Mixins::Loggable::Static::info("Initializing NOX.");
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
      Hermes::Mixins::Loggable::Static::info("Transferring previous fine mesh solution to new fine mesh.");
      Hermes::Hermes2D::OGProjectionNOX<double> ogProjection; ogProjection.project_global(ref_space_new, ref_sln, coeff_vec);
    }

    // Choose preconditioning.
    MlPrecond<double> pc("sa");
    if (PRECOND)
    {
      if (TRILINOS_JFNK) newton_nox.set_precond(pc);
      else newton_nox.set_precond(preconditioner);
    }

    // Perform Newton's iteration.
    Hermes::Mixins::Loggable::Static::info("Performing JFNK on new fine mesh.");
    try
    {
      newton_nox.solve(coeff_vec);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
      
    }

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(newton_nox.get_sln_vector(), ref_space_new, ref_sln);

    // Output.
    Hermes::Mixins::Loggable::Static::info("Number of nonlin iterations: %d (norm of residual: %g)", 
      newton_nox.get_num_iters(), newton_nox.get_residual());
    Hermes::Mixins::Loggable::Static::info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      newton_nox.get_num_lin_iters(), newton_nox.get_achieved_tol());

    // Project the fine mesh solution on the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting fine mesh solution to coarse mesh.");
    Hermes::Hermes2D::OGProjectionNOX<double> ogProjection; ogProjection.project_global(space, ref_sln, sln);

    // Time measurement.
    cpu_time.tick();

    // VTK output.
    if (VTK_VISUALIZATION) 
    {
      // Output solution in VTK format.
      Views::Linearizer lin;
      char* title = new char[100];
      sprintf(title, "sln-%d.vtk", as);
      lin.save_solution_vtk(sln, title, "Potential", false);
      Hermes::Mixins::Loggable::Static::info("Solution in VTK format saved to file %s.", title);

      // Output mesh and element orders in VTK format.
      Views::Orderizer ord;
      sprintf(title, "ord-%d.vtk", as);
      ord.save_orders_vtk(space, title);
      Hermes::Mixins::Loggable::Static::info("Element orders in VTK format saved to file %s.", title);
    }

    // View the coarse mesh solution and polynomial orders.
    if (HERMES_VISUALIZATION) 
    {
      sview.show(sln);
      oview.show(space);
    }

    // Skip visualization time.
    cpu_time.tick();

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      space.get_num_dofs(), ref_space_new->get_num_dofs(), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    cpu_time.tick();    
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    graph_dof.add_values(space.get_num_dofs(), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    
    // Skip the time spent to save the convergence graphs.
    cpu_time.tick();

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) 
      done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(&selector);

      // Increase the counter of performed adaptivity steps.
      if (done == false)  
        as++;
    }

    // Clean up.
    delete [] coeff_vec;
  }
  while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the last fine mesh solution - final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(ref_sln);

  // Wait for all views to be closed.
  Views::View::wait();

  return 0;
}

