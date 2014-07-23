#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// This example shows how to use the Hcurl space and
// visualize finite element basis functions. Note that
// higher-order basis functions in this space comprise
// edge functions associated with mesh edges (tangential
// component is zero on the boundary of the element patch
// associated with the edge), and bubble functions
// associated with elements (tangential component is
// zero on the element boundary).
//
// The following parameters can be changed:

// Initial uniform mesh refinement.
int INIT_REF_NUM = 2;
// Polynomial degree of mesh elements.
int P_INIT = 3;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Create an Hcurl space with default shapeset.
  SpaceSharedPtr<double> space(new HcurlSpace<double>(mesh, P_INIT));

  // Visualize FE basis.
  VectorBaseView<double> bview("VectorBaseView", new WinGeom(0, 0, 700, 600));
  bview.show(space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}