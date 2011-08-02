#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This example is the same as the previous one except that the 
//  nonlinearity is given via a cubic spline.
//
//  PDE: Stationary heat transfer equation with nonlinear thermal
//       conductivity, - div[lambda(u) grad u] + src(x, y) = 0.
//
//  Nonlinearity: cubic spline approximating 1 + pow(u, alpha).
//
//  Domain: square (-10, 10)^2.
//
//  BC: Nonconstant Dirichlet.
//
//  The following parameters can be changed:

const int P_INIT = 2;                             // Initial polynomial degree.
const double NEWTON_TOL = 1e-8;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 4;                   // Number of initial refinements towards boundary.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double heat_src = 1.0;
double alpha = 4.0;

int main(int argc, char* argv[])
{
  // Define nonlinear thermal conductivity lambda(u) via a cubic spline.
  // Step 1: Fill the x values and use lambda_macro(u) = 1 + u^4 for the y values.
#define lambda_macro(x) (1 + Hermes::pow(x, 4))
  Hermes::vector<double> lambda_pts(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0);
  Hermes::vector<double> lambda_val;
  for (unsigned int i = 0; i < lambda_pts.size(); i++) lambda_val.push_back(lambda_macro(lambda_pts[i]));
  // Step 2: Create the cubic spline (and plot it for visual control). 
  double bc_left = 0.0;
  double bc_right = 0.0;
  bool first_der_left = false;
  bool first_der_right = false;
  bool extrapolate_der_left = true;
  bool extrapolate_der_right = true;
  CubicSpline lambda(lambda_pts, lambda_val, bc_left, bc_right, first_der_left, first_der_right,
                     extrapolate_der_left, extrapolate_der_right);
  info("Saving cubic spline into a Pylab file spline.dat.");
  double interval_extension = 3.0; // The interval of definition of the spline will be 
                                   // extended by "interval_extension" on both sides.
  lambda.plot("spline.dat", interval_extension);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential("Bdy");
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof: %d", ndof);

  // Initialize the weak formulation.
  Hermes2DFunction<double> src(-heat_src);
  DefaultWeakFormPoisson<double> wf(HERMES_ANY, &lambda, &src);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Project the initial condition on the FE space to obtain initial 
  // coefficient vector for the Newton's method.
  // NOTE: If you want to start from the zero vector, just define 
  // coeff_vec to be a vector of ndof zeros (no projection is needed).
  info("Projecting to obtain initial vector for the Newton's method.");
  double* coeff_vec = new double[ndof];
  CustomInitialCondition init_sln(&mesh);
  OGProjection<double>::project_global(&space, &init_sln, coeff_vec, matrix_solver); 

  // Initialize Newton solver.
  NewtonSolver<double> newton(&dp, matrix_solver);

  // Perform Newton's iteration.
  bool freeze_jacobian = false;
  if (!newton.solve(coeff_vec, NEWTON_TOL, NEWTON_MAX_ITER, freeze_jacobian)) 
    error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into a Solution.
  Solution<double> sln;
  Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);

  // Get info about time spent during assembling in its respective parts.
  //dp.get_all_profiling_output(std::cout);

  // Clean up.
  delete [] coeff_vec;

  // Visualise the solution and mesh.
  ScalarView<double> s_view("Solution", new WinGeom(0, 0, 440, 350));
  s_view.show_mesh(false);
  s_view.show(&sln);
  OrderView<double> o_view("Mesh", new WinGeom(450, 0, 400, 350));
  o_view.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

