#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;
using namespace Views;

// This example shows how to adapt the mesh to match an arbitrary 
// given function (no PDE solved). 
//
// The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 2;                             
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
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
// H2D_HP_ANISO_P, H2D_HP_ANISO. 
const CandList CAND_LIST = H2D_HP_ANISO;          
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// Stopping criterion for adaptivity.
const double ERR_STOP = 1.0;                      
// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0. 
const double CONV_EXP = 1.0;                              
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);
  
  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, P_INIT);

  // Initialize the weak formulation.
  WeakForm<double> wf_dummy;

  // Initialize coarse and reference mesh solution.
  Solution<double> sln;
  ExactSolutionCustom* ref_sln = NULL;

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

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
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space<double>* ref_space = Space<double>::construct_refined_space(&space);

    // Assign the function f() to the fine mesh.
    info("Assigning f() to the fine mesh.");
    if(ref_sln != NULL) delete ref_sln;
    ref_sln = new ExactSolutionCustom(ref_space->get_mesh());

    // Time measurement.
    cpu_time.tick();
    
    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection<double>::project_global(&space, ref_sln, &sln, matrix_solver); 
   
    // View the coarse mesh solution and polynomial orders.
    sview.show(&sln);
    oview.show(&space);

    // Calculate element errors and total error estimate.
    info("Calculating exact error."); 
    Adapt<double>* adaptivity = new Adapt<double>(&space);
    // Note: the error estimate is now equal to the exact error.
    double err_exact_rel = adaptivity->calc_err_est(&sln, ref_sln) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_exact_rel: %g%%", 
      Space<double>::get_num_dofs(&space), Space<double>::get_num_dofs(ref_space), err_exact_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space<double>::get_num_dofs(&space), err_exact_rel);
    graph_dof.save("conv_dof.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_exact_rel);
    graph_cpu.save("conv_cpu.dat");

    // If err_exact_rel too large, adapt the mesh.
    if (err_exact_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    if (Space<double>::get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete adaptivity;
    if (done == false)
      delete ref_space->get_mesh();
    delete ref_space;
  }
  while (done == false);
  
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show(ref_sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

