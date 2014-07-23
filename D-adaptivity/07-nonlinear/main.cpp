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
// Parameter influencing the candidate selection.
const double THRESHOLD = 0.2;
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1.0;
// Stopping criterion for the Newton's method on coarse mesh.
const double NEWTON_TOL_COARSE = 1e-4;
// Stopping criterion for the Newton's method on fine mesh.
const double NEWTON_TOL_FINE = 1e-8;
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;
// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Selector.
H1ProjBasedSelector<double> selector(CAND_LIST);

// Problem parameters.
double heat_src = 1.0;

int main(int argc, char* argv[])
{
  // Define nonlinear thermal conductivity lambda(u) via a cubic spline.
  // Step 1: Fill the x values and use lambda_macro(u) = 1 + u^4 for the y values.
#define lambda_macro(x) (1 + std::pow(x, 4))
  std::vector<double> lambda_pts({ -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0 });
  std::vector<double> lambda_val;
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
  Hermes::Mixins::Loggable::Static::info("Saving cubic spline into a Pylab file spline.dat.");
  // The interval of definition of the spline will be
  // extended by "interval_extension" on both sides.
  double interval_extension = 3.0;
  lambda.calculate_coeffs();
  lambda.plot("spline.dat", interval_extension);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh->refine_all_elements();
  mesh->refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential("Bdy");
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Set the space to adaptivity.
  adaptivity.set_space(space);

  // Initialize the weak formulation
  Hermes2DFunction<double> f(-heat_src);
  WeakFormSharedPtr<double> wf(new WeakFormsH1::DefaultWeakFormPoisson<double>(HERMES_ANY, &lambda, &f));

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector<double> selector(CAND_LIST, H2DRS_DEFAULT_ORDER);

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>), ref_sln(new Solution<double>);

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  sview.show_mesh(false);
  OrderView oview("Mesh", new WinGeom(450, 0, 400, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  MeshFunctionSharedPtr<double>  init_sln(new InitialSolutionHeatTransfer(mesh));

  // Initialize Newton solver on coarse mesh.
  Hermes::Mixins::Loggable::Static::info("Solving on coarse mesh:");
  NewtonSolver<double> newton_coarse(wf, space);

  try
  {
    newton_coarse.solve(init_sln);
  }
  catch (std::exception& e)
  {
    std::cout << e.what();
  }

  // Translate the resulting coefficient vector into the Solution<double> sln.
  Solution<double>::vector_to_solution(newton_coarse.get_sln_vector(), space, sln);

  // Cleanup after the Newton loop on the coarse mesh.
  NewtonSolver<double> newton(wf, space);
  newton.set_verbose_output(Hermes::Mixins::Loggable::Static::info);

  // Adaptivity loop.
  int as = 1; bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();

    // Initialize discrete problem on the reference mesh.

    // Calculate initial coefficient vector on the reference mesh.
    double* coeff_vec = new double[ref_space->get_num_dofs()];
    if (as == 1)
    {
      // In the first step, project the coarse mesh solution.
      Hermes::Mixins::Loggable::Static::info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
      OGProjection<double>::project_global(ref_space, sln, coeff_vec);
    }
    else
    {
      // In all other steps, project the previous fine mesh solution.
      Hermes::Mixins::Loggable::Static::info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
      OGProjection<double>::project_global(ref_space, ref_sln, coeff_vec);
    }

    // Initialize Newton solver on fine mesh.
    Hermes::Mixins::Loggable::Static::info("Solving on fine mesh:");
    // Perform Newton's iteration.
    newton.set_space(ref_space);
    try
    {
      newton.solve(coeff_vec);
    }
    catch (std::exception& e)
    {
      std::cout << e.what();
    }

    // Translate the resulting coefficient vector into the Solution<double> ref_sln.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution on the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on new coarse mesh for error calculation.");
    OGProjection<double>::project_global(space, ref_sln, sln);

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      Space<double>::get_num_dofs(space), Space<double>::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs(space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");

    // View the coarse mesh solution.
    sview.show(sln);
    oview.show(space);

    // If err_est_rel too large, adapt the mesh.
    if (err_est_rel < ERR_STOP)
      done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting the coarse mesh.");
      done = adaptivity.adapt(&selector);
    }

    // Clean up.
    delete[] coeff_vec;

    as++;
  } while (!done);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(ref_sln);

  // Wait for keyboard or mouse input.
  View::wait();
  return 0;
}