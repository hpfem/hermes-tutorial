#define HERMES_REPORT_ALL
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
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Problem parameters.
const double const_f = 1.0;  

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load(mesh_file, &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_essential("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = Space<double>::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakFormsH1::DefaultWeakFormPoisson<double> wf(HERMES_ANY, 
      new Hermes1DFunction<double>(1.0), new Hermes2DFunction<double>(-const_f));

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
  Vector<double>* rhs = create_vector<double>(matrix_solver);
  LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver, matrix, rhs);

  // Initial coefficient vector for the Newton's method.  
  double* coeff_vec = new double[ndof];
  memset(coeff_vec, 0, ndof*sizeof(double));

  // Perform Newton's iteration.
  Hermes::Hermes2D::Solution<double> sln;
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver);
  try
  {
    newton.solve(coeff_vec);
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.printMsg();
    error("Newton's iteration failed.");
  }

  Hermes::Hermes2D::Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);
  
  // Visualize the solution.
  Views::ScalarView view("Solution", new Views::WinGeom(0, 0, 800, 350));
  view.show(&sln);

  // Wait for the view to be closed.
  Views::View::wait();

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  return 0;
}

