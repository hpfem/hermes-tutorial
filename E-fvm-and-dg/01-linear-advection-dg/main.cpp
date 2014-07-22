#include "definitions.h"

//  This example solves a linear advection equation using Dicontinuous Galerkin (DG) method and standard continuous Finite Element Method. Shows the comparison.
//
//  PDE: \nabla \cdot (\Beta u) = 0, where \Beta = (-x_2, x_1) / |x| represents a circular counterclockwise flow field.
//
//  Domain: Square (0, 1) x (0, 1).
//
//  BC:	Dirichlet, u = 1 where \Beta(x) \cdot n(x) < 0, that is on [0,0.5] x {0}, and g = 0 anywhere else.
//				
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF = 3;
// Number of initial mesh refinements according to a specific criterion.
const int INIT_REF_CRITERION = 4;
// Distance from the arc x^2 + y^2 = 0.25 where to refine.
const double INIT_REF_DIST = 0.07;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 1;
// Use shock capturing for DG.
const bool DG_SHOCK_CAPTURING = true;
// What parts to use.
const bool WANT_DG = true;
const bool WANT_FEM = true;

// Boundary markers.
const std::string BDY_BOTTOM_LEFT = "1";

int criterion(Element * e)
{
  if(std::abs(e->vn[0]->x * e->vn[0]->x + e->vn[0]->y * e->vn[0]->y - 0.25) < INIT_REF_DIST)
    return 0;
  if(std::abs(e->vn[2]->x * e->vn[2]->x + e->vn[2]->y * e->vn[2]->y - 0.25) < INIT_REF_DIST)
    return 0;
  return -1;
}

// Limiting for DG.
#include "euler_util.cpp"

int main(int argc, char* args[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF; i++) 
    mesh->refine_all_elements();

  mesh->refine_by_criterion(criterion, INIT_REF_CRITERION);

  ScalarView view1("Solution - Discontinuous Galerkin FEM", new WinGeom(900, 0, 450, 350));
  ScalarView view2("Solution - Standard continuous FEM", new WinGeom(900, 400, 450, 350));

  if(WANT_DG)
  {
    // Create an L2 space.
    SpaceSharedPtr<double> space_l2(new L2Space<double>(mesh, P_INIT));

    // Initialize the solution.
    MeshFunctionSharedPtr<double> sln_l2(new Solution<double>);

    // Initialize the weak formulation.
    WeakFormSharedPtr<double> wf_l2(new CustomWeakForm(BDY_BOTTOM_LEFT));

    // Initialize linear solver.
    Hermes::Hermes2D::LinearSolver<double> linear_solver(wf_l2, space_l2);

    Hermes::Mixins::Loggable::Static::info("Assembling Discontinuous Galerkin (nelem: %d, ndof: %d).", mesh->get_num_active_elements(), space_l2->get_num_dofs());

    // Solve the linear system. If successful, obtain the solution.
    Hermes::Mixins::Loggable::Static::info("Solving Discontinuous Galerkin.");
    try
    {
      linear_solver.solve();
      if(DG_SHOCK_CAPTURING)
      {      
        FluxLimiter flux_limiter(FluxLimiter::Kuzmin, linear_solver.get_sln_vector(), space_l2, true);

        flux_limiter.limit_second_orders_according_to_detector();

        flux_limiter.limit_according_to_detector();

        flux_limiter.get_limited_solution(sln_l2);

        view1.set_title("Solution - limited Discontinuous Galerkin FEM");
      }
      else
        Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), space_l2, sln_l2);

      // View the solution.
      view1.show(sln_l2);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
}
  }
  if(WANT_FEM)
  {
    // Create an H1 space.
    SpaceSharedPtr<double> space_h1(new H1Space<double>(mesh, P_INIT));

    // Initialize the solution.
    MeshFunctionSharedPtr<double> sln_h1(new Solution<double>);

    // Initialize the weak formulation.
    WeakFormSharedPtr<double> wf_h1(new CustomWeakForm(BDY_BOTTOM_LEFT, false));

    Hermes::Mixins::Loggable::Static::info("Assembling Continuous FEM (nelem: %d, ndof: %d).", mesh->get_num_active_elements(), space_h1->get_num_dofs());

    // Initialize linear solver.
    Hermes::Hermes2D::LinearSolver<double> linear_solver(wf_h1, space_h1);

    // Solve the linear system. If successful, obtain the solution.
    Hermes::Mixins::Loggable::Static::info("Solving Continuous FEM.");
    try
    {
      linear_solver.solve();
      Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), space_h1, sln_h1);

      // View the solution.
      view2.show(sln_h1);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }
  }

  // Wait for keyboard or mouse input.
  View::wait();
  return 0;
}
