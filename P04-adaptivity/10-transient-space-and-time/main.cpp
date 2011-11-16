#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;
using namespace Views;

//  This example is analogous to the previous two, and it combines spatial adaptivity 
//  with adaptive time stepping. If an embedded method is used, temporal error is 
//  measured and visualized. Adaptive time stepping can be turned on or 
//  off using the flag ADAPTIVE_TIME_STEP_ON. An embedded R-K method must be 
//  used if ADAPTIVE_TIME_STEP_ON == true.
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
const int INIT_GLOB_REF_NUM = 2;                   
// Number of initial refinements towards boundary.
const int INIT_BDY_REF_NUM = 0;                    
// Initial polynomial degree.
const int P_INIT = 2;                              
// Time step. 
double time_step = 0.05;                           
// Time interval length.
const double T_FINAL = 5.0;                        

// Spatial adaptivity.
// Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_FREQ = 1;                         
// 1... mesh reset to basemesh and poly degrees to P_INIT.   
// 2... one ref. layer shaved off, poly degrees reset to P_INIT.
// 3... one ref. layer shaved off, poly degrees decreased by one. 
const int UNREF_METHOD = 3;                       
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;                     
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
// More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int STRATEGY = 0;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;          
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0. 
const double CONV_EXP = 1.0;                      
// Stopping criterion for adaptivity.
const double SPACE_ERR_TOL = 1.0;                 
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Matrix solvers: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Temporal adaptivity.
// This flag decides whether adaptive time stepping will be done.
// The methods for the adaptive and fixed-step versions are set
// below. An embedded method must be used with adaptive time stepping. 
bool ADAPTIVE_TIME_STEP_ON = true;                
// If rel. temporal error is greater than this threshold, decrease time 
// step size and repeat time step.
const double TIME_ERR_TOL_UPPER = 1.0;            
// If rel. temporal error is less than this threshold, increase time step
// but do not repeat time step (this might need further research).
const double TIME_ERR_TOL_LOWER = 0.1;            
// Time step increase ratio (applied when rel. temporal error is too small).
const double TIME_STEP_INC_RATIO = 1.1;           
// Time step decrease ratio (applied when rel. temporal error is too large).
const double TIME_STEP_DEC_RATIO = 0.8;           

// Newton's method.
// Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL_COARSE = 0.001;           
// Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL_FINE = 0.005;             
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 20;                   

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
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Turn off adaptive time stepping if R-K method is not embedded.
  if (bt.is_embedded() == false && ADAPTIVE_TIME_STEP_ON == true) {
    warn("R-K method not embedded, turning off adaptive time stepping.");
    ADAPTIVE_TIME_STEP_ON = false;
  }

  // Load the mesh.
  Mesh mesh, basemesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &basemesh);
  mesh.copy(&basemesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  EssentialBCNonConst bc_essential("Bdy");
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();

  // Convert initial condition into a Solution.
  CustomInitialCondition sln_time_prev(&mesh);

  // Initialize the weak formulation
  CustomNonlinearity lambda(alpha);
  Hermes2DFunction<double> f(heat_src);
  WeakFormsH1::DefaultWeakFormPoisson<double> wf(HERMES_ANY, &lambda, &f);

  // Initialize the discrete problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Create a refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Visualize initial condition.
  char title[100];
  ScalarView sln_view("Initial condition", new WinGeom(0, 0, 440, 350));
  sln_view.show_mesh(false);
  OrderView ordview("Initial mesh", new WinGeom(445, 0, 440, 350));
  ScalarView time_error_view("Temporal error", new WinGeom(0, 400, 440, 350));
  time_error_view.fix_scale_width(60);
  ScalarView space_error_view("Spatial error", new WinGeom(445, 400, 440, 350));
  space_error_view.fix_scale_width(60);
  sln_view.show(&sln_time_prev);
  ordview.show(&space);

  // Graph for time step history.
  SimpleGraph time_step_graph;
  if (ADAPTIVE_TIME_STEP_ON) info("Time step history will be saved to file time_step_history.dat.");
  
  // Time stepping loop.
  double current_time = 0.0; int ts = 1;
  do 
  {
    info("Begin time step %d.", ts);
    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
        case 1: mesh.copy(&basemesh);
                space.set_uniform_order(P_INIT);
                break;
        case 2: mesh.unrefine_all_elements();
                space.set_uniform_order(P_INIT);
                break;
        case 3: mesh.unrefine_all_elements();
                space.adjust_element_order(-1, -1, P_INIT, P_INIT);
                break;
        default: error("Wrong global derefinement method.");
      }

      ndof = Space<double>::get_num_dofs(&space);
    }
    info("ndof: %d", ndof);

    // Spatial adaptivity loop. Note: sln_time_prev must not be 
    // changed during spatial adaptivity. 
    Solution<double> ref_sln;
    Solution<double>* time_error_fn;
    if (bt.is_embedded() == true) time_error_fn = new Solution<double>(&mesh);
    else time_error_fn = NULL;
    bool done = false; int as = 1;
    double err_est;
    do {
      // Construct globally refined reference mesh and setup reference space.
      Space<double>* ref_space = Space<double>::construct_refined_space(&space);

      RungeKutta<double> runge_kutta(&wf, ref_space, &bt, matrix_solver);

      // Runge-Kutta step on the fine mesh.
      info("Runge-Kutta time step on fine mesh (t = %g s, tau = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
      bool verbose = true;
      
      try
      {
        runge_kutta.rk_time_step_newton(current_time, time_step, &sln_time_prev, 
                                      &ref_sln, time_error_fn, false, false, verbose, 
                                      NEWTON_TOL_FINE, NEWTON_MAX_ITER);
      }
      catch(Exceptions::Exception& e)
      {
        e.printMsg();
        error("Runge-Kutta time step failed");
      }

      /* If ADAPTIVE_TIME_STEP_ON == true, estimate temporal error. 
         If too large or too small, then adjust it and restart the time step. */

      double rel_err_time = 0;
      if (bt.is_embedded() == true) {
        info("Calculating temporal error estimate.");

        // Show temporal error.
        char title[100];
        sprintf(title, "Temporal error est, spatial adaptivity step %d", as);     
        time_error_view.set_title(title);
        time_error_view.show_mesh(false);
        AbsFilter abs_tef(time_error_fn);
        time_error_view.show(&abs_tef, HERMES_EPS_HIGH);

        rel_err_time = Global<double>::calc_norm(time_error_fn, HERMES_H1_NORM) / 
                       Global<double>::calc_norm(&ref_sln, HERMES_H1_NORM) * 100;
        if (ADAPTIVE_TIME_STEP_ON == false) info("rel_err_time: %g%%", rel_err_time);
      }

      if (ADAPTIVE_TIME_STEP_ON) {
        if (rel_err_time > TIME_ERR_TOL_UPPER) {
          info("rel_err_time %g%% is above upper limit %g%%", rel_err_time, TIME_ERR_TOL_UPPER);
          info("Decreasing tau from %g to %g s and restarting time step.", 
               time_step, time_step * TIME_STEP_DEC_RATIO);
          time_step *= TIME_STEP_DEC_RATIO;
          delete ref_space;
          continue;
        }
        else if (rel_err_time < TIME_ERR_TOL_LOWER) {
          info("rel_err_time = %g%% is below lower limit %g%%", rel_err_time, TIME_ERR_TOL_LOWER);
          info("Increasing tau from %g to %g s.", time_step, time_step * TIME_STEP_INC_RATIO);
          time_step *= TIME_STEP_INC_RATIO;
          delete ref_space;
          continue;
        }
        else {
          info("rel_err_time = %g%% is in acceptable interval (%g%%, %g%%)", 
            rel_err_time, TIME_ERR_TOL_LOWER, TIME_ERR_TOL_UPPER);
        }

        // Add entry to time step history graph.
        time_step_graph.add_values(current_time, time_step);
        time_step_graph.save("time_step_history.dat");
      }

      /* Estimate spatial errors and perform mesh refinement */

      info("Spatial adaptivity step %d.", as);

      // Project the fine mesh solution onto the coarse mesh.
      Solution<double> sln;
      info("Projecting fine mesh solution on coarse mesh for error estimation.");
      OGProjection<double>::project_global(&space, &ref_sln, &sln, matrix_solver); 

      // Show spatial error.
      sprintf(title, "Spatial error est, spatial adaptivity step %d", as);  
      DiffFilter<double>* space_error_fn = new DiffFilter<double>(Hermes::vector<MeshFunction<double>*>(&ref_sln, &sln));   
      space_error_view.set_title(title);
      space_error_view.show_mesh(false);
      AbsFilter abs_sef(space_error_fn);
      space_error_view.show(&abs_sef);

      // Calculate element errors and spatial error estimate.
      info("Calculating spatial error estimate.");
      Adapt<double>* adaptivity = new Adapt<double>(&space);
      double err_rel_space = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

      // Report results.
      info("ndof: %d, ref_ndof: %d, err_rel_space: %g%%", 
           Space<double>::get_num_dofs(&space), Space<double>::get_num_dofs(ref_space), err_rel_space);

      // If err_est too large, adapt the mesh.
      if (err_rel_space < SPACE_ERR_TOL) done = true;
      else 
      {
        info("Adapting the coarse mesh.");
        done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

        if (Space<double>::get_num_dofs(&space) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }
      
      // Clean up.
      delete adaptivity;
      delete ref_space;
      delete space_error_fn;
    }
    while (done == false);

    // Clean up.
    if (time_error_fn != NULL) delete time_error_fn;

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution<double>, time %g s", current_time);
    sln_view.set_title(title);
    sln_view.show_mesh(false);
    sln_view.show(&ref_sln);
    sprintf(title, "Mesh, time %g s", current_time);
    ordview.set_title(title);
    ordview.show(&space);

    // Copy last reference solution into sln_time_prev.
    sln_time_prev.copy(&ref_sln);

    // Increase current time and counter of time steps.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
