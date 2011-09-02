#define HERMES_REPORT_ALL
#include "definitions.h"

// This example shows how to solve exisymmetric problems. The domain of interest
// is a hollow cylinder whose axis is aligned with the y-axis. It has fixed
// temperature on the bottom face, and a radiation (Newton) condition on
// all other faces.
//
// PDE: Laplace equation -div (LAMBDA grad u) = 0 (y-axis is the axis of symmetry),
//      LAMBDA... constant.
//
// BC: u = T1 ... fixed temperature on Gamma_bottom,
//     -LAMBDA * du/dn = ALPHA*(u - T0) ... heat flux on the rest of the boundary (Gamma_heat_flux).
//
// Note that the last BC can be written in the form  du/dn - H*u = -H*T0.
//
// The following parameters can be changed:

const bool USE_XML_FORMAT = true;                 // Select whether you want to read 
                                                  // the original or XML mesh file.
const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int P_INIT = 4;                             // Uniform polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
// Problem parameters.
const double T1 = 100.0;      // Prescribed temperature on Gamma_bottom.
const double T0 = 0.0;        // Outer temperature.
const double LAMBDA = 386;    // Thermal conductivity.
const double ALPHA = 20.0;    // Heat_flux coefficient on Gamma_heat_flux.

// Weak forms.
#include "definitions.h"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  if (USE_XML_FORMAT == true)
  {
    MeshReaderH2DXML mloader;  
    info("Reading mesh in XML format.");
    mloader.load("domain.xml", &mesh);
  }
  else 
  {
    MeshReaderH2D mloader;
    info("Reading mesh in original format.");
    mloader.load("domain.mesh", &mesh);
  }

  // Perform initial mesh refinements.
  for(int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential("Bottom", T1);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakFormPoissonNewton wf(LAMBDA, ALPHA, T0, "Heat_flux");

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Initial coefficient vector for the Newton's method.  
  double* coeff_vec = new double[ndof];
  memset(coeff_vec, 0, ndof*sizeof(double));

  // Initialize Newton solver.
  NewtonSolver<double> newton(&dp, matrix_solver);

  // Perform Newton's iteration.
  if (!newton.solve(coeff_vec)) 
    error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into a Solution.
  Solution<double> sln;
  Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);

  // VTK output.
  if (VTK_VISUALIZATION) 
  {
    // Output solution in VTK format.
    Linearizer lin(&sln);
    bool mode_3D = true;
    lin.save_solution_vtk("sln.vtk", "Temperature", mode_3D);
    info("Solution in VTK format saved to file %s.", "sln.vtk");

    // Output mesh and element orders in VTK format.
    Orderizer ord;
    ord.save_orders_vtk(&space, "ord.vtk");
    info("Element orders in VTK format saved to file %s.", "ord.vtk");
  }

  // Visualize the solution.
  if (HERMES_VISUALIZATION) 
  {
    ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
    // Hermes uses adaptive FEM to approximate higher-order FE solutions with linear
    // triangles for OpenGL. The second parameter of View::show() sets the error 
    // tolerance for that. Options are HERMES_EPS_LOW, HERMES_EPS_NORMAL (default), 
    // HERMES_EPS_HIGH and HERMES_EPS_VERYHIGH. The size of the graphics file grows 
    // considerably with more accurate representation, so use it wisely.
    view.show(&sln, HERMES_EPS_HIGH);
    View::wait();
  }

  // Clean up.
  delete [] coeff_vec;

  return 0;
}

