#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;
using namespace Views;

//  This example shows how the Newton's method can be combined with
//  automatic adaptivity.
//
//  PDE: stationary heat transfer equation with nonlinear thermal
//  conductivity, -div[S(u)grad u] - heat_src = 0 where S(u) is a cubic spline.
//
//  Domain: unit square (-10,10)^2.
//
//  BC: Dirichlet.
//
//  The following parameters can be changed:

// Initial polynomial degree.
const int P_INIT = 1;                             
// Number of initial uniform mesh refinements.
const int INIT_GLOB_REF_NUM = 1;                  
// Number of initial refinements towards boundary.
const int INIT_BDY_REF_NUM = 3;                   
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.2;                     
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
// More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int STRATEGY = 1;                           
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
const double ERR_STOP = 1.0;                      
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Stopping criterion for the Newton's method on coarse mesh.
const double NEWTON_TOL_COARSE = 1e-4;            
// Stopping criterion for the Newton's method on fine mesh.
const double NEWTON_TOL_FINE = 1e-8;              
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                  
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Problem parameters.
double heat_src = 1.0;

int main(int argc, char* argv[])
{

  // Define nonlinear thermal conductivity lambda(u) via a cubic spline.
  // Step 1: Fill the x values and use lambda_macro(u) = 1 + u^4 for the y values.
  #define lambda_macro(x) (1 + std::pow(x, 4))
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
  // The interval of definition of the spline will be 
  // extended by "interval_extension" on both sides.
  double interval_extension = 3.0; 
  lambda.plot("spline.dat", interval_extension);

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

  // Initialize the weak formulation
  Hermes2DFunction<double> f(-heat_src);
  WeakFormsH1::DefaultWeakFormPoisson<double> wf(HERMES_ANY, &lambda, &f);

  // Initialize the FE problem.
  DiscreteProblem<double> dp_coarse(&wf, &space);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize coarse and reference mesh solution.
  Solution<double> sln, ref_sln;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  sview.show_mesh(false);
  OrderView oview("Mesh", new WinGeom(450, 0, 400, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector on the coarse mesh.");
  double* coeff_vec_coarse = new double[space.get_num_dofs()];
  InitialSolutionHeatTransfer init_sln(&mesh);
  OGProjection<double>::project_global(&space, &init_sln, coeff_vec_coarse, matrix_solver);

  // Initialize Newton solver on coarse mesh.
  info("Solving on coarse mesh:");
  bool verbose = true;
  NewtonSolver<double> newton_coarse(&dp_coarse, matrix_solver);
  newton_coarse.set_verbose_output(verbose);

  // Perform initial Newton's iteration on coarse mesh, to obtain 
  // good initial guess for the Newton's method on the fine mesh.
  try
  {
    newton_coarse.solve(coeff_vec_coarse, NEWTON_TOL_COARSE, NEWTON_MAX_ITER);
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.printMsg();
    error("Newton's iteration failed.");
  }

  // Translate the resulting coefficient vector into the Solution<double> sln.
  Solution<double>::vector_to_solution(newton_coarse.get_sln_vector(), &space, &sln);

  // Cleanup after the Newton loop on the coarse mesh.
  delete [] coeff_vec_coarse;

  // Adaptivity loop.
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space<double>* ref_space = Space<double>::construct_refined_space(&space);

    // Initialize discrete problem on the reference mesh.
    DiscreteProblem<double> dp(&wf, ref_space);

    // Calculate initial coefficient vector on the reference mesh.
    double* coeff_vec = new double[ref_space->get_num_dofs()];
    if (as == 1)
    {
      // In the first step, project the coarse mesh solution.
      info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
      OGProjection<double>::project_global(ref_space, &sln, coeff_vec, matrix_solver);
    }
    else
    {
      // In all other steps, project the previous fine mesh solution.
      info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
      OGProjection<double>::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver);
      delete ref_sln.get_space();
      delete ref_sln.get_mesh();
    }

    // Initialize Newton solver on fine mesh.
    info("Solving on fine mesh:");
    bool verbose = true;
    NewtonSolver<double> newton(&dp, matrix_solver);
    newton.set_verbose_output(verbose);

    // Perform Newton's iteration.
    try
    {
      newton.solve(coeff_vec, NEWTON_TOL_FINE, NEWTON_MAX_ITER);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    }

    // Translate the resulting coefficient vector into the Solution<double> ref_sln.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);

    // Project the fine mesh solution on the coarse mesh.
    info("Projecting reference solution on new coarse mesh for error calculation.");
    OGProjection<double>::project_global(&space, &ref_sln, &sln, matrix_solver);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate.");
    Adapt<double>* adaptivity = new Adapt<double>(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      Space<double>::get_num_dofs(&space), Space<double>::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs(&space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");

    // View the coarse mesh solution.
    sview.show(&sln);
    oview.show(&space);

    // If err_est_rel too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      info("Adapting the coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (Space<double>::get_num_dofs(&space) >= NDOF_STOP)
      {
        done = true;
        break;
      }
    }

    // Clean up.
    delete [] coeff_vec;
    delete adaptivity;

    as++;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(&ref_sln);

  // Wait for keyboard or mouse input.
  View::wait();
  return 0;
}

