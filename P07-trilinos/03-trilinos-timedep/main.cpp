#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

//  The purpose of this example is to show how to use Trilinos
//  for time-dependent PDE problem.
//  NOX solver is used, either using Newton's method or JFNK and
//  with or without preconditioning,
//
//  PDE: Heat transfer: HEATCAP*RHO*du/dt - div(LAMBDA * grad u) = 0.
//
//  Domain: Unit square.
//
//  BC: Dirichlet at the bottom, Newton du/dn = ALPHA*(TEMP_EXT - u) elsewhere.
//

const int INIT_REF_NUM = 4;       // Number of initial uniform mesh refinements.
const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const double ALPHA = 10.0;        // Coefficient for the Nwwton boundary condition.
const double LAMBDA = 1e5;
const double HEATCAP = 1e6;
const double RHO = 3000.0;
const double TEMP_EXT = 20.0;
const double TEMP_INIT = 10.0;

const double TAU = 50.0;          // Time step.        

const bool JFNK = true;
const bool PRECOND = true;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc("Bdy_bottom", TEMP_INIT);
  EssentialBCs<double> bcs(&bc);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = Space<double>::get_num_dofs(&space);
  info("ndof: %d", ndof);

  // Define constant initial condition. 
  ConstantSolution<double> t_prev_time(&mesh, TEMP_INIT);

  // Initialize the weak formulation.
  CustomWeakForm wf(Hermes::vector<std::string>("Bdy_right", "Bdy_top", "Bdy_left"), HEATCAP, RHO, TAU, LAMBDA, ALPHA, TEMP_EXT, &t_prev_time, JFNK);

  // Initialize the finite element problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Project the function "t_prev_time" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution on the FE mesh.");
  double* coeff_vec = new double[ndof];
  OGProjection<double>::project_global(&space, &t_prev_time, coeff_vec);

  // Initialize NOX solver.
  NoxSolver<double> solver(&dp);

  // Select preconditioner.
  MlPrecond<double> pc("sa");
  if (PRECOND)
  {
    if (JFNK) solver.set_precond(pc);
    else solver.set_precond("ML");
  }

  // Initialize the view.
  ScalarView Tview("Temperature", new WinGeom(0, 0, 450, 400));
  Tview.set_min_max_range(10,20);

  // Time stepping loop:
  double total_time = 0.0;
  for (int ts = 1; total_time <= 2000.0; ts++)
  {
    info("---- Time step %d, t = %g s", ts, total_time += TAU);

    info("Assembling by DiscreteProblem, solving by NOX.");
    if (solver.solve(coeff_vec))
      Solution<double>::vector_to_solution(solver.get_sln_vector(), &space, &t_prev_time);
    else
      error("NOX failed.");

    // Show the new solution.
    Tview.show(&t_prev_time);

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      solver.get_num_lin_iters(), solver.get_achieved_tol());
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
