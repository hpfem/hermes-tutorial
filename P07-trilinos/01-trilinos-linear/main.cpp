#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  The purpose of this example is to use a simple linear problem with known 
//  exact solution to show how to use NOX, and to compare its performance to 
//  other linear solvers in Hermes (MUMPS, PETSc, UMFPACK, etc.).
//
//  PDE: Poisson equation.
//
//  Domain: Square (-1, 1)^2.
//
//  BC: Nonhomogeneous Dirichlet, see the function essential_bc_values() below.
//
//  Exact solution: x*x + y*y.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 6;                       // Number of initial uniform mesh refinements.
const int P_INIT = 3;                             // Initial polynomial degree of all mesh elements.
const bool TRILINOS_JFNK = true;                  // true = Jacobian-free method (for NOX),
                                                  // false = Newton (for NOX).
const bool PRECOND = true;                        // Preconditioning by jacobian in case of JFNK (for NOX),
                                                  // default ML preconditioner in case of Newton.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "GMRES";              // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "AztecOO";            // Name of the preconditioner employed by AztecOO 
                                                  // Possibilities: None" - No preconditioning. 
                                                  // "AztecOO" - AztecOO internal preconditioner.
                                                  // "New Ifpack" - Ifpack internal preconditioner.
                                                  // "ML" - Multi level preconditione
// NOX parameters.
unsigned message_type = NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration | NOX::Utils::InnerIteration | NOX::Utils::Parameters | NOX::Utils::LinearSolverDetails;
                                                  // NOX error messages, see NOX_Utils.h.

double ls_tolerance = 1e-5;                       // Tolerance for linear system.
unsigned flag_absresid = 0;                       // Flag for absolute value of the residuum.
double abs_resid = 1.0e-3;                        // Tolerance for absolute value of the residuum.
unsigned flag_relresid = 1;                       // Flag for relative value of the residuum.
double rel_resid = 1.0e-2;                        // Tolerance for relative value of the residuum.
int max_iters = 100;                              // Max number of iterations.

int main(int argc, char **argv)
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Set exact solution.
  CustomExactSolution exact(&mesh);

  // Initialize boundary conditions
  DefaultEssentialBCNonConst<double> bc_essential("Bdy", &exact);
  EssentialBCs<double> bcs(&bc_essential);
  
  // Initialize the weak formulation.
  CustomWeakFormPoisson wf1;
 
  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = Space<double>::get_num_dofs(&space);
  info("ndof: %d", ndof);

  info("---- Assembling by DiscreteProblem, solving by %s:", 
       MatrixSolverNames[matrix_solver_type].c_str());

  // Initialize the linear discrete problem.
  DiscreteProblem<double> dp1(&wf1, &space);
    
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver_type);
  Vector<double>* rhs = create_vector<double>(matrix_solver_type);
  LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver_type, matrix, rhs);
  
  #ifdef HAVE_AZTECOO
    if (matrix_solver_type == SOLVER_AZTECOO) 
    {
      (dynamic_cast<AztecOOSolver<double>*>(solver))->set_solver(iterative_method);
      (dynamic_cast<AztecOOSolver<double>*>(solver))->set_precond(preconditioner);
      // Using default iteration parameters (see solver/aztecoo.h).
    }    
#endif
  
  // Begin time measurement of assembly.
  cpu_time.tick(HERMES_SKIP);

  // Initial coefficient vector for the Newton's method.  
  double* coeff_vec = new double[ndof];
  memset(coeff_vec, 0, ndof*sizeof(double));

  // Initialize the Newton solver.
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp1, matrix_solver_type);

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::Solution<double> sln1;
  if (!newton.solve(coeff_vec)) 
    error("Newton's iteration failed.");
  else
    Hermes::Hermes2D::Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln1);

  // CPU time measurement.
  double time = cpu_time.tick().last();

  // Calculate errors.
  double rel_err_1 = Global<double>::calc_rel_error(&sln1, &exact, HERMES_H1_NORM) * 100;
  info("CPU time: %g s.", time);
  info("Exact H1 error: %g%%.", rel_err_1);

  delete(matrix);
  delete(rhs);
  delete(solver);
    
  // View the solution and mesh.
  ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  sview.show(&sln1);
  //OrderView  oview("Polynomial orders", new WinGeom(450, 0, 400, 350));
  //oview.show(&space);
  
  // TRILINOS PART:

  info("---- Assembling by DiscreteProblem, solving by NOX:");

  // Initialize the weak formulation for Trilinos.
  CustomWeakFormPoisson wf2(TRILINOS_JFNK);
  
  // Initialize DiscreteProblem.
  DiscreteProblem<double> dp2(&wf2, &space);
  
  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Set initial vector for NOX.
  // NOTE: Using zero vector was causing convergence problems.
  info("Projecting to obtain initial vector for the Newton's method.");
  ZeroSolution init_sln(&mesh);

  // Initialize the NOX solver with the vector "coeff_vec".
  info("Initializing NOX.");
  NoxSolver<double> nox_solver(&dp2);
  nox_solver.set_output_flags(message_type);

  nox_solver.set_ls_type(iterative_method);
  nox_solver.set_ls_tolerance(ls_tolerance);

  nox_solver.set_conv_iters(max_iters);
  if (flag_absresid)
    nox_solver.set_conv_abs_resid(abs_resid);
  if (flag_relresid)
    nox_solver.set_conv_rel_resid(rel_resid);

  // Choose preconditioning.
  MlPrecond<double> pc("sa");
  if (PRECOND)
  {
    if (TRILINOS_JFNK) nox_solver.set_precond(pc);
    else nox_solver.set_precond(preconditioner);
  }

  // Assemble and solve using NOX.
  Solution<double> sln2;
  OGProjection<double>::project_global(&space, &init_sln, coeff_vec);
  if (nox_solver.solve(coeff_vec))
  {
    Solution<double>::vector_to_solution(nox_solver.get_sln_vector(), &space, &sln2);

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      nox_solver.get_num_iters(), nox_solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      nox_solver.get_num_lin_iters(), nox_solver.get_achieved_tol());
  }
  else error("NOX failed");

  // CPU time needed by NOX.
  time = cpu_time.tick().last();

  // Show the NOX solution.
  ScalarView view2("Solution<double> 2", new WinGeom(450, 0, 460, 350));
  view2.show(&sln2);
  //view2.show(&exact);

  // Calculate errors.
  double rel_err_2 = Global<double>::calc_rel_error(&sln2, &exact, HERMES_H1_NORM) * 100;
  info("CPU time: %g s.", time);
  info("Exact H1 error: %g%%.", rel_err_2);
 
  // Wait for all views to be closed.
  View::wait();
  
  return 0;
}
