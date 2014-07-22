#include "definitions.h"

//  This example solves a general second-order linear equation with non-constant
//  coefficients, and shows how integration orders in linear and bilinear forms
//  can be defined manually.
//
//  PDE: -d/dx(a_11(x,y)du/dx) - d/dx(a_12(x,y)du/dy) - d/dy(a_21(x,y)du/dx) - d/dy(a_22(x,y)du/dy)
//       + a_1(x,y)du/dx + a_2(x,y)du/dy + a_0(x,y)u - rhs(x,y) = 0.
//
//  Domain: arbitrary
//
//  BC:  Dirichlet for boundary marker "Horizontal": u = g_D(x,y)
//       Natural for boundary marker "Vertical":   (a_11(x,y)*nu_1 + a_21(x,y)*nu_2) * dudx
//                                               + (a_12(x,y)*nu_1 + a_22(x,y)*nu_2) * dudy = g_N(x,y)
//
//  The following parameters can be changed:

// Read the original or XML mesh file.
const bool USE_XML_FORMAT = true;                 
// Initial polynomial degree of all mesh elements.
const int P_INIT = 3;                             
// Number of initial uniform refinements.
const int INIT_REF_NUM = 3;

int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

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

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();
  
  // Initialize boundary conditions
  CustomEssentialBCNonConst bc_essential("Horizontal");
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();
  Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf(new CustomWeakFormGeneral("Horizontal"));

  // Initialize the FE problem.
  DiscreteProblem<double> dp(wf, space);

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

  // Time measurement.
  cpu_time.tick();

  // View the solution and mesh.
  ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  sview.show(sln);
  OrderView oview("Polynomial orders", new WinGeom(450, 0, 405, 350));
  oview.show(space);

  // Print timing information.
  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();

  return 0;
}
