#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// This example illustrates how to use full-featured NURBS. Simplified
// format is enabled for circular arcs (see example 03-poisson). 
//
// PDE: Poisson equation -Laplace u - const_f = 0 with homogeneous (zero)
//      Dirichlet boundary conditions.
//
// Domain: Rectangle (0, 2) x (0, 1) where the upper edge is a NURBS
//         (see the end of the mesh file for details).

// Choose one of the following mesh files:
// Mesh containing NURBS with one control point.
//const char* mesh_file = "domain-1.mesh";            
// Mesh containing NURBS with two control points.
const char* mesh_file = "domain-2.mesh";          
// Mesh containing NURBS with three control points.
//const char* mesh_file = "domain-3.mesh";          

// The following parameters can be also changed:
// Uniform polynomial degree of mesh elements.
const int P_INIT = 3;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;

// Problem parameters.
const double const_f = 1.0;  

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load(mesh_file, mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_essential("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = Space<double>::get_num_dofs(space);
  Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakFormsH1::DefaultWeakFormPoisson<double> wf(HERMES_ANY, 
      new Hermes1DFunction<double>(1.0), new Hermes2DFunction<double>(-const_f));

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, space);

  // Perform Newton's iteration.
  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  NewtonSolver<double> newton(&dp);
  try
  {
    newton.solve();
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
    
  }

  Hermes::Hermes2D::Solution<double>::vector_to_solution(newton.get_sln_vector(), space, sln);
  
  // Visualize the solution.
  Views::ScalarView view("Solution", new Views::WinGeom(0, 0, 800, 350));
  view.show(sln);

  // Wait for the view to be closed.
  Views::View::wait();

  return 0;
}

