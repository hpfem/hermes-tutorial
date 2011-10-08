#define HERMES_REPORT_ALL
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
const bool WANT_FEM = false;
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;

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
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF; i++) 
    mesh.refine_all_elements();

  mesh.refine_by_criterion(criterion, INIT_REF_CRITERION);

  MeshView m;
  m.show(&mesh);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver_type);
  Vector<double>* rhs = create_vector<double>(matrix_solver_type);
  LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver_type, matrix, rhs);

  ScalarView view1("Solution - Discontinuous Galerkin FEM", new WinGeom(900, 0, 450, 350));
  ScalarView view2("Solution - Standard continuous FEM", new WinGeom(900, 400, 450, 350));

  if(WANT_DG)
  {
    // Create an L2 space.
    L2Space<double> space_l2(&mesh, P_INIT);
  
    // Initialize the solution.
    Solution<double> sln_l2;

    // Initialize the weak formulation.
    CustomWeakForm wf_l2(BDY_BOTTOM_LEFT);


    // Initialize the FE problem.
    DiscreteProblem<double> dp_l2(&wf_l2, &space_l2);
    
    info("Assembling Discontinuous Galerkin (nelem: %d, ndof: %d).", mesh.get_num_active_elements(), space_l2.get_num_dofs());
    dp_l2.assemble(matrix, rhs);

    // Solve the linear system. If successful, obtain the solution.
    info("Solving Discontinuous Galerkin.");
    if(solver->solve())
      if(DG_SHOCK_CAPTURING)
      {      
        FluxLimiter flux_limiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), &space_l2, true);

        flux_limiter.limit_second_orders_according_to_detector();

        flux_limiter.limit_according_to_detector();

        flux_limiter.get_limited_solution(&sln_l2);

        view1.set_title("Solution - limited Discontinuous Galerkin FEM");
      }
      else
        Solution<double>::vector_to_solution(solver->get_sln_vector(), &space_l2, &sln_l2);
    else 
      error ("Matrix solver failed.\n");

    // View the solution.
    view1.show(&sln_l2);
  }
  if(WANT_FEM)
  {
    // Create an H1 space.
    H1Space<double> space_h1(&mesh, P_INIT);
  
    // Initialize the solution.
    Solution<double> sln_h1;

    // Initialize the weak formulation.
    CustomWeakForm wf_h1(BDY_BOTTOM_LEFT, false);


    // Initialize the FE problem.
    DiscreteProblem<double> dp_h1(&wf_h1, &space_h1);
    
    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver_type);
    Vector<double>* rhs = create_vector<double>(matrix_solver_type);
    LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver_type, matrix, rhs);

    info("Assembling Continuous FEM (nelem: %d, ndof: %d).", mesh.get_num_active_elements(), space_h1.get_num_dofs());
    dp_h1.assemble(matrix, rhs);

    // Solve the linear system. If successful, obtain the solution.
    info("Solving Continuous FEM.");
    if(solver->solve())
      Solution<double>::vector_to_solution(solver->get_sln_vector(), &space_h1, &sln_h1);
    else 
      error ("Matrix solver failed.\n");
    
    // View the solution.
    view2.show(&sln_h1);
  }

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;
  
  // Wait for keyboard or mouse input.
  View::wait();
  return 0;
}
