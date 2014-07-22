#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

//  The purpose of this example is to show how to use Trilinos for nonlinear PDE problems. It 
//  compares performance of the Newton's method in Hermes (assembling and solving via the NewtonSolver<double> 
//  class and matrix problem solution via UMFpack) with the performance of the Trilinos/NOX 
//  solver (using the internal Hermes DiscreteProblem<double> class to assemble discrete problems).
//
//  PDE:  - \nabla (k \nabla u) - f = 0
//  k = (1 + sqr(u_x) + sqr(u_y))^{-0.5}
//
//  Domain: Unit square.
//
//  BC: zero Dirichlet.
//
//  Exact solution: (x - x*x) * (y - y*y).
//
//  Initial guess for the Newton's method: zero function.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 5;                       
// Initial polynomial degree of all mesh elements.
const int P_INIT = 3;                             
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-6;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                  
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// NOX parameters.
// true = jacobian-free method,
// false = Newton.
const bool TRILINOS_JFNK = false;                 
// Preconditioning by jacobian (1) or approximation of jacobian (2)
// in case of JFNK,
// Default ML proconditioner in case of Newton.
const int PRECOND = 2;                            
// Name of the iterative method employed by AztecOO (ignored
// by the other solvers). 
// Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* iterative_method = "bicgstab";        
// Name of the preconditioner employed by AztecOO (ignored by
// the other solvers).
// Possibilities: none, jacobi, neumann, least-squares, or a
//  preconditioner from IFPACK (see solver/aztecoo.h)
const char* preconditioner = "least-squares";     
// NOX error messages, see NOX_Utils.h.
unsigned message_type = NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration | NOX::Utils::InnerIteration | NOX::Utils::Parameters | NOX::Utils::LinearSolverDetails;
                                                  
// Tolerance for linear system.
double ls_tolerance = 1e-5;                       
// Flag for absolute value of the residuum.
unsigned flag_absresid = 0;                       
// Tolerance for absolute value of the residuum.
double abs_resid = 1.0e-8;                        
// Flag for relative value of the residuum.
unsigned flag_relresid = 1;                       
// Tolerance for relative value of the residuum.
double rel_resid = 1.0e-8;                        
// Max number of iterations.
int max_iters = 100;                              

int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = Space<double>::get_num_dofs(space);
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

  Hermes::Mixins::Loggable::Static::info("Assembling by DiscreteProblem, solving by Umfpack:");

  // Time measurement.
  cpu_time.tick();

  // Initialize weak formulation,
  WeakFormSharedPtr<double> wf1(new CustomWeakForm);

  // Initialize the Newton solver.
  Hermes::Hermes2D::NewtonSolver<double> newton(wf1, space);

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  MeshFunctionSharedPtr<double> sln1(new Solution<double>);
  MeshFunctionSharedPtr<double> sln2(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_tmp(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double> exact(new CustomExactSolution(mesh));

  try
  {
    newton.solve();
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
}

  // Translate the solution vector into a Solution.
  Hermes::Hermes2D::Solution<double>::vector_to_solution(newton.get_sln_vector(), space, sln1);

  // CPU time needed by UMFpack
  double time1 = cpu_time.tick().last();

  // Time measurement.
  cpu_time.tick();
 
  // Show UMFPACK solution.
  ScalarView view1("Solution 1", new WinGeom(0, 0, 500, 400));
  view1.show(sln1);

  // Calculate error.
  DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator1(RelativeErrorToGlobalNorm, 1);
  errorCalculator1.calculate_errors(sln1, exact);
  double rel_err_1 = errorCalculator1.get_total_error_squared() * 100;

  Hermes::Mixins::Loggable::Static::info("Solution 1 :  exact H1 error: %g%% (time %g s)", rel_err_1, time1);

  // TRILINOS PART:

  // Project the initial condition to obtain the initial
  // coefficient vector.
  Hermes::Mixins::Loggable::Static::info("Projecting to obtain initial vector for the Newton's method.");
  double* coeff_vec = new double[ndof];
  OGProjection<double>::project_global(space, sln_tmp, coeff_vec);

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize the weak formulation for Trilinos.
  WeakFormSharedPtr<double> wf2(new CustomWeakForm(TRILINOS_JFNK, PRECOND == 1, PRECOND == 2));

  // Initialize the NOX solver with the vector "coeff_vec".
  Hermes::Mixins::Loggable::Static::info("Initializing NOX.");
  NewtonSolverNOX<double> solver_nox(wf2, space);
  solver_nox.set_output_flags(message_type);

  solver_nox.set_ls_tolerance(ls_tolerance);

  solver_nox.set_conv_iters(max_iters);
  if (flag_absresid)
    solver_nox.set_conv_abs_resid(abs_resid);
  if (flag_relresid)
    solver_nox.set_conv_rel_resid(rel_resid);

  // Choose preconditioning.
  MlPrecond<double> pc("sa");
  if (PRECOND)
  {
    if (TRILINOS_JFNK) solver_nox.set_precond(pc);
    else solver_nox.set_precond("ML");
  }

  // Solve the nonlinear problem using NOX.
  Hermes::Mixins::Loggable::Static::info("Assembling by DiscreteProblem, solving by NOX.");
  try
  {
    solver_nox.solve(coeff_vec);
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
}

  Solution<double>::vector_to_solution(solver_nox.get_sln_vector(), space, sln2);
  Hermes::Mixins::Loggable::Static::info("Number of nonlin iterations: %d (norm of residual: %g)", 
        solver_nox.get_num_iters(), solver_nox.get_residual());
  Hermes::Mixins::Loggable::Static::info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
        solver_nox.get_num_lin_iters(), solver_nox.get_achieved_tol());

  // CPU time needed by NOX.
  double time2 = cpu_time.tick().last();

  // Calculate error.
  DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator2(RelativeErrorToGlobalNorm, 1);
  errorCalculator2.calculate_errors(sln2, exact);
  double rel_err_2 = errorCalculator2.get_total_error_squared() * 100;
  Hermes::Mixins::Loggable::Static::info("Solution 2 (NOX): exact H1 error: %g%% (time %g + %g = %g [s])", rel_err_2, proj_time, time2, proj_time+time2);

  // Show NOX solution.
  ScalarView view2("Solution 2", new WinGeom(510, 0, 500, 400));
  view2.show(sln2);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
