#include "definitions.h"

using namespace Hermes::Hermes2D::Views;

// This example shows how to use the L2 finite element space and L2 shapeset.
// As a sample problem, a continuous function x^3 + y^3 is projected onto the
// L2 finite element space in the L2 norm. When zero-order is used, the result
// is a piecewice constant function. The class BaseView will show you the basis
// functions.
//
// The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 4;
// Polynomial degree of mesh elements.
const int P_INIT = 1;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform uniform mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Create an L2 space with default shapeset.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, P_INIT));

  // View basis functions.
  BaseView<double> bview("BaseView", new WinGeom(0, 0, 600, 500));
  bview.show(space);
  // View::wait(H2DV_WAIT_KEYPRESS);

  // Initialize the exact and projected solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_exact(new CustomExactSolution(mesh));

  // Project the exact function on the FE space.
  OGProjection<double> ogProjection;
  ogProjection.project_global(space, sln_exact, sln);

  // Visualize the projection.
  ScalarView view1("Projection", new WinGeom(610, 0, 600, 500));
  view1.show(sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}