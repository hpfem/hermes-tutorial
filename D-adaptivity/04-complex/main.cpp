#include "definitions.h"

using namespace Hermes::Hermes2D::RefinementSelectors;

typedef std::complex<double> complex;

//  This problem describes the distribution of the vector potential in
//  a 2D domain comprising a wire carrying electrical current, air, and
//  an iron which is not under voltage.
//
//  PDE: -(1/mu)Laplace A + ii*omega*gamma*A - J_ext = 0.
//
//  Domain: Rectangle of height 0.003 and width 0.004. Different
//  materials for the wire, air, and iron (see mesh file domain2.mesh).
//
//  BC: Zero Dirichlet on the top and right edges, zero Neumann
//  elsewhere.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Initial polynomial degree of mesh elements.
const int P_INIT = 1;                             
// Parameter influencing the candidate selection.

const double THRESHOLD = 0.3;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;          
// Stopping criterion for adaptivity.
const double ERR_STOP = 1.0;                      
// Error calculation & adaptivity.
DefaultErrorCalculator<::complex, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<::complex> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<::complex> adaptivity(&errorCalculator, &stoppingCriterion);
// Selector.
H1ProjBasedSelector<::complex> selector(CAND_LIST);

// Problem parameters.
const double MU_0 = 4.0*M_PI*1e-7;
const double MU_IRON = 1e3 * MU_0;
const double GAMMA_IRON = 6e6;
const double J_EXT = 1e6;
const double FREQ = 5e3;
const double OMEGA = 2 * M_PI * FREQ;

int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Initialize boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<::complex> 
      bc_essential("Dirichlet", ::complex(0.0, 0.0));
  EssentialBCs<::complex> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<::complex> space(new  H1Space<::complex>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();
  Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakForm wf("Air", MU_0, "Iron", MU_IRON, GAMMA_IRON,
    "Wire", MU_0, ::complex(J_EXT, 0.0), OMEGA);

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<::complex> sln(new Solution<::complex>), ref_sln(new Solution<::complex>);

  // Initialize refinement selector.
  H1ProjBasedSelector<::complex> 
      selector(CAND_LIST, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 600, 350));
  sview.show_mesh(false);
  Views::OrderView oview("Polynomial orders", new Views::WinGeom(610, 0, 520, 350));

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<::complex>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<::complex> ref_space = ref_space_creator.create_ref_space();
    int ndof_ref = ref_space->get_num_dofs();

    // Initialize reference problem.
    Hermes::Mixins::Loggable::Static::info("Solving on reference mesh.");
    DiscreteProblem<::complex> dp(&wf, ref_space);

    // Time measurement.
    cpu_time.tick();

    // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
    Hermes::Hermes2D::NewtonSolver<::complex> newton(&dp);

    try{
      newton.solve();
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
      
    }
    Hermes::Hermes2D::Solution<::complex>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Time measurement.
    cpu_time.tick();

    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
    OGProjection<::complex> ogProjection; ogProjection.project_global(space, ref_sln, sln);

    // View the coarse mesh solution and polynomial orders.
    MeshFunctionSharedPtr<double> real_filter(new RealFilter(sln));
    sview.show(real_filter);

    oview.show(space);

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      space->get_num_dofs(), ref_space->get_num_dofs(), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(space->get_num_dofs(), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(&selector);
    }

    // Increase counter.
    as++;
  }
  while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");

  RealFilter real_filter(ref_sln);
  sview.show(&real_filter);

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
