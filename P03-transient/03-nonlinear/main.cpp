#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;

//  This example solves a nonlinear time-dependent problem using arbitrary 
//  Runge-Kutta methods.
//
//  For a list of available R-K methods see the file hermes_common/tables.h.
//
//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity, du/dt = div[lambda(u) grad u] + f.
//
//  Nonlinearity: lambda(u) = 1 + pow(u, alpha).
//
//  Domain: square (-10, 10)^2.
//
//  BC: Nonconstant Dirichlet.
//
//  IC: Custom initial condition matching the BC.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_GLOB_REF_NUM = 3;                   
// Number of initial refinements towards boundary.
const int INIT_BDY_REF_NUM = 4;                    
// Initial polynomial degree.
const int P_INIT = 2;                              
// Time step.
const double time_step = 0.2;                      
// Time interval length.
const double T_FINAL = 5.0;                        
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-5;                    
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                   
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   

// Choose one of the following time-integration methods, or define your own Butcher's table. The last number
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods:
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2,
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4,
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_5_4.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
//   Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded,
//   Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded,
//   Implicit_DIRK_ISMAIL_7_45_embedded.
ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

// Problem parameters.
// Parameter for nonlinear thermal conductivity.
const double alpha = 4.0;                         
const double heat_src = 1.0;

// Main function.
int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  EssentialBCNonConst bc_essential("Bdy");
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d.", ndof);

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition sln_time_prev(&mesh);

  // Initialize the weak formulation
  CustomNonlinearity lambda(alpha);
  Hermes2DFunction<double> f(heat_src);
  DefaultWeakFormPoisson<double> wf(HERMES_ANY, &lambda, &f);

  // Next time level solution.
  Solution<double> sln_time_new(&mesh);

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(&wf, &space, &bt, matrix_solver);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 500, 400));
  OrderView oview("Mesh", new WinGeom(510, 0, 460, 400));
  oview.show(&space);

  // Time stepping loop:
  double current_time = 0; int ts = 1;
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).",
         current_time, time_step, bt.get_size());
    bool freeze_jacobian = false;
    bool block_diagonal_jacobian = false;
    bool verbose = true;
    double damping_coeff = 1.0;
    double max_allowed_residual_norm = 1e10;
    Hermes::vector<Solution<double>*> slns_time_prev;
    slns_time_prev.push_back(&sln_time_prev);
    Hermes::vector<Solution<double>*> slns_time_new;
    slns_time_new.push_back(&sln_time_new);
    try
    {
      runge_kutta.rk_time_step_newton(current_time, time_step, slns_time_prev, slns_time_new, 
          freeze_jacobian, block_diagonal_jacobian, verbose, NEWTON_TOL, 
          NEWTON_MAX_ITER, damping_coeff, max_allowed_residual_norm);
    }
    catch(Exceptions::Exception& e)
    {
      e.printMsg();
      error("Runge-Kutta time step failed");
    }

    // Update time.
    current_time += time_step;

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Solution, t = %g", current_time);
    sview.set_title(title);
    sview.show(&sln_time_new, HERMES_EPS_VERYHIGH);
    oview.show(&space);

    // Copy solution for the new time step.
    sln_time_prev.copy(&sln_time_new);

    // Increase counter of time steps.
    ts++;
  }
  while (current_time < T_FINAL);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
