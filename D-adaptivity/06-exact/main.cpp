#include "definitions.h"

using namespace RefinementSelectors;
using namespace Views;

// This example shows how to adapt the mesh to match an arbitrary 
// given function (no PDE solved). 
//
// The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 2;                             
// Parameter influencing the candidate selection.
const double THRESHOLD = 0.2;                        
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
// H2D_HP_ANISO_P, H2D_HP_ANISO. 
const CandList CAND_LIST = H2D_HP_ANISO;  
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-2;
// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Selector.
H1ProjBasedSelector<double> selector(CAND_LIST);

int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);
  
  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, P_INIT));
  
  // Set the space to adaptivity.
  adaptivity.set_space(space);

  // Initialize the weak formulation.
  WeakForm<double> wf_dummy;

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  MeshFunctionSharedPtr<double> ref_sln(new ExactSolutionCustom(mesh));
   
  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView sview("Scalar potential Phi", new WinGeom(0, 0, 610, 300));
  sview.fix_scale_width(40);
  sview.show_mesh(false);
  OrderView  oview("Mesh", new WinGeom(620, 0, 600, 300));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

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

    // Assign the function f() to the fine mesh.
    Hermes::Mixins::Loggable::Static::info("Assigning f() to the fine mesh.");

    // Time measurement.
    cpu_time.tick();
    
    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
    OGProjection<double>::project_global(space, ref_sln, sln); 
   
    // View the coarse mesh solution and polynomial orders.
    sview.show(sln);
    oview.show(space);

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating exact error."); 
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_exact_rel = errorCalculator.get_total_error_squared() * 100;

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_exact_rel: %g%%", 
      Space<double>::get_num_dofs(space), Space<double>::get_num_dofs(ref_space), err_exact_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space<double>::get_num_dofs(space), err_exact_rel);
    graph_dof.save("conv_dof.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_exact_rel);
    graph_cpu.save("conv_cpu.dat");

    // If err_exact_rel too large, adapt the mesh.
    if (err_exact_rel < ERR_STOP) done = true;
    else 
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(&selector);
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
  }
  while (done == false);
  
  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show(ref_sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}