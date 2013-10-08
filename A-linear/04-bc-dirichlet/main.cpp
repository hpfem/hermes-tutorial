#include "definitions.h"

// This example shows how to use non-constant Dirichlet boundary conditions.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = A*x + B*y + C where the 
//                      constants A, B, C are arbitrary constants (defined below).
//
// The following parameters can be changed:

// Read original or XML mesh file.
const bool USE_XML_FORMAT = true;                 
// Set to "false" to suppress Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;           
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;              
// Uniform polynomial degree of mesh elements.
const int P_INIT = 5;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;

// Problem parameters.
// Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_AL = 236.0;            
// Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            
// Volume heat sources generated by electric current.        
const double VOLUME_HEAT_SRC = 6e3;        
const double BDY_A_PARAM = 1.0;
const double BDY_B_PARAM = 2.0;
const double BDY_C_PARAM = 20.0;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  if (USE_XML_FORMAT == true)
  {
    MeshReaderH2DXML mloader;  
    Hermes::Mixins::Loggable::Static::info("Reading mesh in XML format.");
    mloader.load("domain.xml", mesh);
  }
  else 
  {
    MeshReaderH2D mloader;
    Hermes::Mixins::Loggable::Static::info("Reading mesh in original format.");
    mloader.load("domain.mesh", mesh);
  }

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakFormPoissonDirichlet wf("Aluminum", LAMBDA_AL, "Copper", 
                                    LAMBDA_CU, VOLUME_HEAT_SRC);
  
  // Initialize boundary conditions.
  CustomDirichletCondition bc_essential(
      Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"),
      BDY_A_PARAM, BDY_B_PARAM, BDY_C_PARAM);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();
  Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, space);

  // Initialize Newton solver.
  NewtonSolver<double> newton(&dp);

  // Perform Newton's iteration.
  try
  {
    newton.solve();
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
    
  }

  // Translate the resulting coefficient vector into a Solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  Solution<double>::vector_to_solution(newton.get_sln_vector(), space, sln);

  // VTK output.
  if (VTK_VISUALIZATION) 
  {
    // Output solution in VTK format.
    Linearizer lin;
    bool mode_3D = true;
    lin.save_solution_vtk(sln, "sln.vtk", "Temperature", mode_3D);
    Hermes::Mixins::Loggable::Static::info("Solution in VTK format saved to file %s.", "sln.vtk");

    // Output mesh and element orders in VTK format.
    Orderizer ord;
    ord.save_orders_vtk(space, "ord.vtk");
    Hermes::Mixins::Loggable::Static::info("Element orders in VTK format saved to file %s.", "ord.vtk");
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
    view.show(sln, HERMES_EPS_HIGH);
    View::wait();
  }

  return 0;
}

