#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This example uses the Picard's method to solve a nonlinear problem.
//
//  PDE: Stationary heat transfer equation with nonlinear thermal
//  conductivity, -div[lambda(u) grad u] + src(x, y) = 0.
//
//  Nonlinearity: lambda(u) = 1 + Hermes::pow(u, 4).
//
//  Picard's linearization: -div[lambda(u^n) grad u^{n+1}] + src(x, y) = 0.
//
//  Domain: square (-10, 10)^2.
//
//  BC: Nonconstant Dirichlet.
//
//  The following parameters can be changed:

const int P_INIT = 2;                             // Initial polynomial degree.
const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 5;                   // Number of initial refinements towards boundary.
const double PICARD_TOL = 1e-2;                   // Stopping criterion for the Picard's method.
const int PICARD_MAX_ITER = 1000;                 // Maximum allowed number of Picard's iterations.
const double INIT_COND_CONST = 3.0;               // Value for custom constant initial condition.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double heat_src = 1.0;
double alpha = 4.0;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
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

  // Initialize previous iteration solution for the Picard's method.
  CustomInitialCondition init_sln(&mesh, INIT_COND_CONST);

  // Project the initial condition on the FE space to obtain initial 
  // coefficient vector for the Picard's method.
  // NOTE: If you want to start from the zero vector, just define 
  // coeff_vec to be a vector of ndof zeros (no projection is needed).
  info("Projecting to obtain initial vector for the Picard's method.");
  double* coeff_vec = new double[ndof];
  OGProjection<double>::project_global(&space, &init_sln, coeff_vec, matrix_solver); 

  // FIXME - converting the coefficient vector back to a Solution.
  Solution<double> sln_prev_iter;
  Solution<double>::vector_to_solution(coeff_vec, &space, &sln_prev_iter);

  // Initialize the weak formulation.
  CustomNonlinearity lambda(alpha);
  Hermes2DFunction<double> src(-heat_src);
  CustomWeakFormPicard wf(&sln_prev_iter, &lambda, &src);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Initialize the Picard solver.
  PicardSolver<double> picard(&dp, &sln_prev_iter, matrix_solver);

  // Perform the Picard's iteration.
  Solution<double> sln;
  int number_of_last_iterations_used = 4;
  double anderson_beta = 1.0;
    if (!picard.solve(PICARD_TOL, PICARD_MAX_ITER, number_of_last_iterations_used, anderson_beta)) 
    error("Picard's iteration failed.");
  else
    Solution<double>::vector_to_solution(picard.get_sln_vector(), &space, &sln);
  
  // Visualise the solution and mesh.
  ScalarView<double> s_view("Solution", new WinGeom(0, 0, 440, 350));
  s_view.show_mesh(false);
  s_view.show(&sln_prev_iter);
  OrderView<double> o_view("Mesh", new WinGeom(450, 0, 420, 350));
  o_view.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

