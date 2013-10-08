#include "definitions.h"

using namespace RefinementSelectors;

// This example shows how to run adaptive hp-FEM, h-FEM and p-FEM with
// basic control parameters. The underlying problem is a planar model
// of an electrostatic micromotor (MEMS). You may want to experiment with
// various types of adaptivity via the options H2D_P_ISO, H2D_P_ANISO,
// H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H, H2D_HP_ANISO_P,
// and H2D_HP_ANISO. See the User Documentation for more details.
//
// Uniform initial polynomial degree of mesh elements can be set using
// the variable P_INIT. Before using adaptivity, you have to define a refinement
// selector as shown below. The function adapt() takes the selector
// as a parameter, along with THRESHOLD, STRATEGY, and MESH_REGULARITY.
//   
// Additional control parameters are available, these will be demonstrated
// in the following tutorial examples. In this example, two types of convergence
// graphs are created -- error estimate wrt. the number of degrees of freedom
// (DOF), and error estimate wrt. CPU time. Later we will show how to output
// the error wrt. exact solution when exact solution is available.
//   
// This example also demonstrates how to define different material parameters
// in various parts of the computational domain, and how to measure time.
//
// PDE: -div[eps_r(x,y) grad phi] = 0
//      eps_r = EPS_1 in Omega_1 (surrounding air)
//      eps_r = EPS_2 in Omega_2 (moving part of the motor)
//
// BC: phi = 0 V on Gamma_1 (left edge and also the rest of the outer boundary
//     phi = VOLTAGE on Gamma_2 (boundary of stator)
//
// The following parameters can be changed:

// Set to "false" to suppress Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;           
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;             
// Initial polynomial degree of mesh elements.
const int P_INIT = 1;                             
// Parameter influencing the candidate selection.
const double THRESHOLD = 0.8;                     
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
// H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO_H;        
// Stopping criterion for adaptivity.
const double ERR_STOP = 0.5;                      
// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Selector.
H1ProjBasedSelector<double> selector(CAND_LIST);

// Problem parameters.
const double EPS0 = 8.863e-12;
const double VOLTAGE = 50.0;
const double EPS_MOTOR = 10.0 * EPS0;
const double EPS_AIR = 1.0 * EPS0;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Motor", EPS_MOTOR, "Air", EPS_AIR);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential_out("Outer", 0.0);
  DefaultEssentialBCConst<double> bc_essential_stator("Stator", VOLTAGE);
  EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_essential_out, &bc_essential_stator));

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Initialize coarse and fine mesh solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>), ref_sln(new Solution<double>);

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 410, 600));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  Views::OrderView  oview("Polynomial orders", new Views::WinGeom(420, 0, 400, 600));

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;

  DiscreteProblem<double> dp(&wf, space);
  NewtonSolver<double> newton(&dp);
  newton.set_verbose_output(true);

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);
    
    // Time measurement.
    cpu_time.tick();

    // Construct globally refined mesh and setup fine mesh space.
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();
    int ndof_ref = ref_space->get_num_dofs();

    // Initialize fine mesh problem.
    Hermes::Mixins::Loggable::Static::info("Solving on fine mesh.");
    
    newton.set_space(ref_space);


    // Perform Newton's iteration.
    try
    {
      newton.solve();
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);
    
    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting fine mesh solution on coarse mesh.");
    OGProjection<double>::project_global(space, ref_sln, sln);

    // Time measurement.
    cpu_time.tick();

    // VTK output.
    if (VTK_VISUALIZATION) 
    {
      // Output solution in VTK format.
      Views::Linearizer lin;
      char* title = new char[100];
      sprintf(title, "sln-%d.vtk", as);
      lin.save_solution_vtk(sln, title, "Potential", false);
      Hermes::Mixins::Loggable::Static::info("Solution in VTK format saved to file %s.", title);

      // Output mesh and element orders in VTK format.
      Views::Orderizer ord;
      sprintf(title, "ord-%d.vtk", as);
      ord.save_orders_vtk(space, title);
      Hermes::Mixins::Loggable::Static::info("Element orders in VTK format saved to file %s.", title);
    }

    // View the coarse mesh solution and polynomial orders.
    if (HERMES_VISUALIZATION) 
    {
      sview.show(sln);
      oview.show(space);
    }

    // Skip visualization time.
    cpu_time.tick();

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
    Adapt<double> adaptivity(space);
    bool solutions_for_adapt = true;
    // In the following function, the Boolean parameter "solutions_for_adapt" determines whether
    // the calculated errors are intended for use with adaptivity (this may not be the case, for example,
    // when error wrt. an exact solution is calculated). The default value is solutions_for_adapt = true,
    // The last parameter "error_flags" determine whether the total and element errors are treated as
    // absolute or relative. Its default value is error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL.
    // In subsequent examples and benchmarks, these two parameters will be often used with
    // their default values, and thus they will not be present in the code explicitly.
    double err_est_rel = adaptivity.calc_err_est(sln, ref_sln, solutions_for_adapt,
                         HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      space->get_num_dofs(), ref_space->get_num_dofs(), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    cpu_time.tick();    
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    graph_dof.add_values(space->get_num_dofs(), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    
    // Skip the time spent to save the convergence graphs.
    cpu_time.tick();

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) 
      done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of performed adaptivity steps.
      if (done == false)  
        as++;
    }
    if (space->get_num_dofs() >= NDOF_STOP) 
      done = true;
  }
  while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the fine mesh solution - final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(ref_sln);

  // Wait for all views to be closed.
  Views::View::wait();

  return 0;
}

