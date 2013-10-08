#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// This example shows how to load a mesh, perform various types
// of initial refinements, and use keyboard and mouse controls.
//
// Geometry: L-Shape domain (see file domain.mesh and domain.xml).

// Read original or XML mesh file.
const bool USE_XML_FORMAT = true;     

// Text message with hints.
static char text[] = "\
Click into the image window and:\n\
  press 'm' to show/hide element material markers,\n\
  press 'i' to show/hide element indices,\n\
  press 'b' to toggle boundary markers,\n\
  enlarge your window and press 'c' to center the mesh,\n\
  zoom into the mesh using the right mouse button\n\
  move the mesh around using the left mouse button\n\
  press 'c' to center the mesh again,\n\
  press 'h' to render high resolution image,\n\
  press 's' to save a screenshot in bmp format\n\
  press 'q' to quit.\n\
  Press F1 for help.\n";

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  if (USE_XML_FORMAT == true)
  {
    MeshReaderH2DXML mloader;  
    Hermes::Mixins::Loggable::Static::info("Reading mesh in XML format.");
    try
    {
      mloader.load("domain.xml", mesh);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.print_msg();
    }
  }
  else 
  {
    MeshReaderH2D mloader;
    Hermes::Mixins::Loggable::Static::info("Reading mesh in original format.");
    mloader.load("domain.mesh", mesh);
  }

  // Refine mesh uniformly (optional).
  mesh->refine_all_elements();          

  // Refine towards a mesh vertex (optional).
  // Four refinements towards vertex no. 3.  
  mesh->refine_towards_vertex(3, 4);    

  // Refine towards boundary (optional).
  // Four successive refinements towards boundary with marker "Outer".
  mesh->refine_towards_boundary("Outer", 4);  

  // Refine individual elements (optional).
  // 0... isotropic refinement.
  mesh->refine_element_id(86, 0);          
  // 0... isotropic refinement.
  mesh->refine_element_id(112, 0);         
  // 2... anisotropic refinement.
  mesh->refine_element_id(84, 2);          
  // 1... anisotropic refinement.
  mesh->refine_element_id(114, 1);         

  // Display the mesh.
  // (0, 0) is the upper left corner position, 
  // 350 x 350 is the window size.
  MeshView mview("Hello world!", new WinGeom(0, 0, 350, 350));
  mview.show(mesh);

  // Practice some keyboard and mouse controls.
  printf("%s", text);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
