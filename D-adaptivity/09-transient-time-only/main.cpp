#include "definitions.h"

using namespace RefinementSelectors;
using namespace Views;

//  This example is analogous to example 07-transient-space-only except that 
//  a fixed mesh is used and only the time stepping is adaptive. An arbitrary 
//  embedded Runge-Kutta method can be used. By embedded we mean that the 
//  Butcher's table contains two B rows. The two B rows are used to calculate 
//  two different approximations Y_{n+1} on the next time level, with different 
//  orders of accuracy. With those one works as usual.
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
double time_step = 0.5;                           
// Time interval length.
const double T_FINAL = 5.0;                        
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-5;                    
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                   
// If rel. temporal error is greater than this threshold, decrease time 
// step size and repeat time step.
const double TIME_TOL_UPPER = 1.0;                 
// If rel. temporal error is less than this threshold, increase time step
// but do not repeat time step (this might need further research).
const double TIME_TOL_LOWER = 0.5;                 
// Time step increase ratio (applied when rel. temporal error is too small).
const double TIME_STEP_INC_RATIO = 2.0;            
// Time step decrease ratio (applied when rel. temporal error is too large).
const double TIME_STEP_DEC_RATIO = 0.5;  

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
ButcherTableType butcher_table_type = Implicit_SDIRK_CASH_3_23_embedded;

// Problem parameters.
// Parameter for nonlinear thermal conductivity.
const double alpha = 4.0;                         
const double heat_src = 1.0;

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh->refine_all_elements();
  mesh->refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  EssentialBCNonConst bc_essential("Bdy");
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();

  // Convert initial condition into a Solution.
  MeshFunctionSharedPtr<double>  sln_time_prev(new CustomInitialCondition(mesh));
  MeshFunctionSharedPtr<double>  sln_time_new(new Solution<double>(mesh));
  MeshFunctionSharedPtr<double>  time_error_fn(new ZeroSolution<double>(mesh));

  // Initialize the weak formulation
  CustomNonlinearity lambda(alpha);
  Hermes2DFunction<double> f(heat_src);
  WeakFormSharedPtr<double> wf(new WeakFormsH1::DefaultWeakFormPoisson<double>(HERMES_ANY, &lambda, &f));

  // Initialize views.
  ScalarView sview_high("Solution (higher-order)", new WinGeom(0, 0, 500, 400));
  ScalarView eview("Temporal error", new WinGeom(500, 0, 500, 400));
  eview.fix_scale_width(50);

  RungeKutta<double> runge_kutta(wf, space, &bt);

  // Graph for time step history.
  SimpleGraph time_step_graph;
  Hermes::Mixins::Loggable::Static::info("Time step history will be saved to file time_step_history.dat.");

  // Time stepping loop:
  double current_time = 0.0; int ts = 1;
  do 
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    Hermes::Mixins::Loggable::Static::info("Runge-Kutta time step (t = %g, tau = %g, stages: %d).", 
         current_time, time_step, bt.get_size());
    try
    {
      runge_kutta.set_verbose_output(true);
      runge_kutta.set_time(current_time);
      runge_kutta.set_time_step(time_step);
      runge_kutta.set_max_allowed_iterations(NEWTON_MAX_ITER);
      runge_kutta.set_tolerance(NEWTON_TOL);
      runge_kutta.rk_time_step_newton(sln_time_prev, sln_time_new, time_error_fn);
    }
    catch(Exceptions::Exception& e)
    {
      std::cout << e.what();
    }

    // Plot error function.
    char title[100];
    sprintf(title, "Temporal error, t = %g", current_time);
    eview.set_title(title);
    MeshFunctionSharedPtr<double> abs_tef(new AbsFilter(time_error_fn));
    eview.show(abs_tef);
    
    // Show the new time level solution.
    sprintf(title, "Solution (higher-order), t = %g", current_time);
    sview_high.set_title(title);
    sview_high.show(sln_time_new);

    DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculatorTime(RelativeErrorToGlobalNorm, 1);
    errorCalculatorTime.calculate_errors(time_error_fn, sln_time_new);
    double rel_err_time = errorCalculatorTime.get_total_error_squared() * 100;


    Hermes::Mixins::Loggable::Static::info("rel_err_time = %g%%", rel_err_time);
    if (rel_err_time > TIME_TOL_UPPER) {
      Hermes::Mixins::Loggable::Static::info("rel_err_time above upper limit %g%% -> decreasing time step from %g to %g and repeating time step.", 
           TIME_TOL_UPPER, time_step, time_step * TIME_STEP_DEC_RATIO);
      time_step *= TIME_STEP_DEC_RATIO;
      continue;
    }
    if (rel_err_time < TIME_TOL_LOWER) {
      Hermes::Mixins::Loggable::Static::info("rel_err_time = below lower limit %g%% -> increasing time step from %g to %g", 
           TIME_TOL_UPPER, time_step, time_step * TIME_STEP_INC_RATIO);
      time_step *= TIME_STEP_INC_RATIO;
    }
   
    // Add entry to the timestep graph.
    time_step_graph.add_values(current_time, time_step);
    time_step_graph.save("time_step_history.dat");

    // Copy solution for next time step.
    sln_time_prev->copy(sln_time_new);

    // Update time.
    current_time += time_step;

    // Increase counter of time steps.
    ts++;
} 
  while (current_time < T_FINAL);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
