#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  The purpose of this example is to show how to use Trilinos for nonlinear time-dependent coupled PDE systems.
//  Solved by NOX solver via Newton or JFNK, with or without preconditioning.
//
//  PDE: Flame propagation (same as tutorial example 19-newton-timedep-flame).
//
//  Domain: Same as in tutorial example 19-newton-timedep-flame.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;         
// Initial polynomial degree of all mesh elements.
const int P_INIT = 1;               
// Time step.
const double TAU = 0.5;            
// Time interval length.
const double T_FINAL = 60.0;        

// Problem parameters.
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

// NOX parameters.
// Display more details about nonlinear and linear solvers.
const bool TRILINOS_OUTPUT = true;                
// true = Jacobian-free method (for NOX),
// false = Newton (for NOX).
const bool TRILINOS_JFNK = true;
// Preconditioning by jacobian (1) (less GMRES iterations, more time to create precond)
// or by approximation of jacobian (2) (less time for precond creation, more GMRES iters).
// in case of jfnk, default Ifpack proconditioner in case of Newton.
const int PRECOND = 1;                            
// Name of the iterative method employed by AztecOO (ignored by the other solvers). 
// Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* iterative_method = "GMRES";           
// Name of the preconditioner employed by AztecOO 
// Possibilities: None" - No preconditioning. 
// "AztecOO" - AztecOO internal preconditioner.
// "New Ifpack" - Ifpack internal preconditioner.
// "ML" - Multi level preconditione
const char* preconditioner = "ML";           
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
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> left_t("Left", 1.0);
  EssentialBCs<double> bcs_t(&left_t);

  DefaultEssentialBCConst<double> left_c("Left", 0.0);
  EssentialBCs<double> bcs_c(&left_c);

  // Create H1 spaces with default shapesets.
  H1Space<double>* t_space = new H1Space<double>(&mesh, &bcs_t, P_INIT);
  H1Space<double>* c_space = new H1Space<double>(&mesh, &bcs_c, P_INIT);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(t_space, c_space));
  info("ndof = %d.", ndof);

  // Define initial conditions.
  InitialSolutionTemperature t_prev_time_1(&mesh, x1);
  InitialSolutionConcentration c_prev_time_1(&mesh, x1, Le);
  InitialSolutionTemperature t_prev_time_2(&mesh, x1);
  InitialSolutionConcentration c_prev_time_2(&mesh, x1, Le);
  Solution<double> t_prev_newton;
  Solution<double> c_prev_newton;

  // Filters for the reaction rate omega and its derivatives.
  CustomFilter omega(Hermes::vector<Solution<double>*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);
  CustomFilterDt omega_dt(Hermes::vector<Solution<double>*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);
  CustomFilterDc omega_dc(Hermes::vector<Solution<double>*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);

  // Initialize visualization.
  ScalarView rview("Reaction rate", new WinGeom(0, 0, 800, 230));

  // Initialize weak formulation.
  CustomWeakForm wf(Le, alpha, beta, kappa, x1, TAU, TRILINOS_JFNK, PRECOND, &omega, &omega_dt, 
                    &omega_dc, &t_prev_time_1, &c_prev_time_1, &t_prev_time_2, &c_prev_time_2);

  // Project the functions "t_prev_time_1" and "c_prev_time_1" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solutions on the FE meshes.");
  double* coeff_vec = new double[ndof];
  OGProjection<double>::project_global(Hermes::vector<const Space<double> *>(t_space, c_space), 
                                       Hermes::vector<MeshFunction<double>*>(&t_prev_time_1, &c_prev_time_1),
                                       coeff_vec);

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize finite element problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<const Space<double>*>(t_space, c_space));

  // Initialize NOX solver and preconditioner.
  NewtonSolverNOX<double> solver(&dp);
  MlPrecond<double> pc("sa");
  if (PRECOND)
  {
    if (TRILINOS_JFNK) solver.set_precond(pc);
    else solver.set_precond("Ifpack");
  }
  if (TRILINOS_OUTPUT)
    solver.set_output_flags(NOX::Utils::Error | NOX::Utils::OuterIteration |
                            NOX::Utils::OuterIterationStatusTest |
                            NOX::Utils::LinearSolverDetails);

  // Time stepping loop:
  double total_time = 0.0;
  cpu_time.tick_reset();
  for (int ts = 1; total_time <= T_FINAL; ts++)
  {
    info("---- Time step %d, t = %g s", ts, total_time + TAU);

    cpu_time.tick(HERMES_SKIP);
    try
    {
      solver.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("NOX failed.");
    }

    Solution<double>::vector_to_solutions(solver.get_sln_vector(), Hermes::vector<const Space<double> *>(t_space, c_space), 
              Hermes::vector<Solution<double> *>(&t_prev_newton, &c_prev_newton));

    cpu_time.tick();
    info("Number of nonlin iterations: %d (norm of residual: %g)",
        solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)",
        solver.get_num_lin_iters(), solver.get_achieved_tol());

    // Time measurement.
    cpu_time.tick(HERMES_SKIP);
			
    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Update global time.
    total_time += TAU;

    // Saving solutions for the next time step.
    if(ts > 1)
    {
      t_prev_time_2.copy(&t_prev_time_1);
      c_prev_time_2.copy(&c_prev_time_1);
    }

    t_prev_time_1.copy(&t_prev_newton);
    c_prev_time_1.copy(&c_prev_newton);
      
    // Visualization.
    rview.set_min_max_range(0.0,2.0);
    rview.show(&omega);
    cpu_time.tick(HERMES_SKIP);

    info("Total running time for time level %d: %g s.", ts, cpu_time.tick().last());
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

