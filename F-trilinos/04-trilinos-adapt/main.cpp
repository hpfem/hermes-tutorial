#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Preconditioners;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

//  The purpose of this example is to show how to use Trilinos while adapting mesh.
//  Solved by NOX solver, either using Newton's method or JFNK, with or without
//  preconditioning. The underlying problem is benchmark "layer-interior".
//
//  PDE: -Laplace u - f = 0.
//
//  Domain: Unit square.
//
//  BC: Nonhomogeneous Dirichlet.
//
//  Known exact solution, see definitions.cpp.
//
//  The following parameters can be changed:

// Initial polynomial degree of all mesh elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// Parameter influencing the candidate selection.
const double THRESHOLD = 0.3;
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
// See User Documentation for details.
const CandList CAND_LIST = H2D_HP_ANISO_H;
// Stopping criterion for adaptivity (rel. error tolerance between the
// fine mesh and coarse mesh solution in percent).
const double ERR_STOP = 26.0;

// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Selector.
H1ProjBasedSelector<double> selector(CAND_LIST);

// Problem parameters.
// Slope of the layer inside the domain.
double slope = 60;

// NOX parameters.
// true = Jacobian-free method (for NOX),
// false = Newton (for NOX).
const bool TRILINOS_JFNK = true;
// Preconditioning by jacobian in case of JFNK (for NOX),
// default ML preconditioner in case of Newton.
const bool PRECOND = true;
// Name of the iterative method employed by AztecOO (ignored
// by the other solvers).
// Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* iterative_method = "GMRES";
// Name of the preconditioner employed by AztecOO
// Possibilities: None" - No preconditioning.
// "AztecOO" - AztecOO internal preconditioner.
// "New Ifpack" - Ifpack internal preconditioner.
// "ML" - Multi level preconditioner.
const char* preconditioner = "AztecOO";
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
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);     // quadrilaterals

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Define exact solution.
  MeshFunctionSharedPtr<double> exact(new CustomExactSolution(mesh, slope));

  // Define right-hand side.
  CustomFunction f(slope);

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf(new WeakFormsH1::DefaultWeakFormPoisson<double>(HERMES_ANY, new Hermes1DFunction<double>(1.0), &f));

  // Initialize boundary conditions
  DefaultEssentialBCNonConst<double> bc_essential("Bdy", exact);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  sview.show_mesh(false);
  sview.fix_scale_width(50);
  OrderView  oview("Polynomial orders", new WinGeom(450, 0, 420, 350));
  OrderView  oviewa("Polynomial orders", new WinGeom(450, 0, 420, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu, graph_dof_exact, graph_cpu_exact;

  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  MeshFunctionSharedPtr<double> ref_sln(new Solution<double>);

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();
    int ndof_ref = Space<double>::get_num_dofs(ref_space);

    // Initialize reference problem.
    Hermes::Mixins::Loggable::Static::info("Solving on reference mesh.");

    // Time measurement.
    cpu_time.tick();

    // Initial coefficient vector for the Newton's method.
    double* coeff_vec = new double[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(double));

    // Initialize NOX solver.
    NewtonSolverNOX<double> solver(wf, ref_space);
    solver.set_output_flags(message_type);

    solver.set_ls_tolerance(ls_tolerance);

    solver.set_conv_iters(max_iters);
    if (flag_absresid)
      solver.set_conv_abs_resid(abs_resid);
    if (flag_relresid)
      solver.set_conv_rel_resid(rel_resid);

    // Select preconditioner.
    MlPrecond<double> pc("sa");
    if (PRECOND)
    {
      if (TRILINOS_JFNK) solver.set_precond(pc);
      else solver.set_precond("ML");
    }

    // Time measurement.
    cpu_time.tick();

    Hermes::Mixins::Loggable::Static::info("Assembling by DiscreteProblem, solving by NOX.");
    try
    {
      solver.solve(coeff_vec);
    }
    catch (std::exception& e)
    {
      std::cout << e.what();
    }

    Solution<double>::vector_to_solution(solver.get_sln_vector(), ref_space, ref_sln);
    Hermes::Mixins::Loggable::Static::info("Number of nonlin iterations: %d (norm of residual: %g)",
      solver.get_num_iters(), solver.get_residual());
    Hermes::Mixins::Loggable::Static::info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)",
      solver.get_num_lin_iters(), solver.get_achieved_tol());

    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
    OGProjection<double>::project_global(space, ref_sln, sln);

    // View the coarse mesh solution and polynomial orders.
    sview.show(sln);
    oview.show(space);
    oviewa.show(ref_space);

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate and exact error.");
    errorCalculator.calculate_errors(sln, exact, false);
    double err_exact_rel = errorCalculator.get_total_error_squared() * 100;
    // Calculate exact error.
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d",
      Space<double>::get_num_dofs(space), Space<double>::get_num_dofs(ref_space));
    Hermes::Mixins::Loggable::Static::info("err_est_rel: %g%%, err_exact_rel: %g%%", err_est_rel, err_exact_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space<double>::get_num_dofs(space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    graph_dof_exact.add_values(Space<double>::get_num_dofs(space), err_exact_rel);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel);
    graph_cpu_exact.save("conv_cpu_exact.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(&selector);

      // Increase the counter of adaptivity steps.
      if (!done)  as++;
    }

    // Clean up.
    delete[] coeff_vec;
  } while (!done);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}