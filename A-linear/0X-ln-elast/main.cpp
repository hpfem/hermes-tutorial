#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
//#include "function/function.h"
#include <fstream>

// Read original or XML mesh file.
const bool USE_XML_FORMAT = true;
// Initial polynomial degree of all elements.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 5;
//
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;

// Problem parameters.
// Young modulus (soft tissue).
const double E  = 1000;
// Poisson ratio.
const double nu = 0.49;
// Non constant lambda (vol Lame coeff)
CustomLambda lambdaFun(E,nu);
// Non constant mu (shear Lame coeff)
CustomMu muFun(E,nu);
// Density.
const double rho = 8000.0;
// Gravitational acceleration.
const double g1 = -9.81;
// Surface force in x-direction.
const double f0  = 0;
// Surface force in y-direction.
const double f1  = 8e4;


int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, mesh1;
  if (USE_XML_FORMAT == true)
  {
    MeshReaderH2DXML mloader;
    Hermes::Mixins::Loggable::Static::info("Reading mesh in XML format.");
    mloader.load("square.xml", &mesh);
  }
  else
  {
    MeshReaderH2D mloader;
    Hermes::Mixins::Loggable::Static::info("Reading mesh in original format.");
    mloader.load("square.mesh", &mesh);
  }

  // Perform uniform mesh refinement.
  int refinement_type = 0;
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(refinement_type);

   // Show mesh.
   MeshView mv("Mesh", new WinGeom(0, 0, 580, 400));
   mv.show(&mesh);

   // Exact lambda
   ExactSolutionLambda exact_lambda(&mesh,E,nu);
   ScalarView viewLam("lambda [Pa]", new WinGeom(0, 460, 530, 350));
   viewLam.show_mesh(false);
   viewLam.show(&exact_lambda);
   // Exact lambda
   ExactSolutionMu exact_mu(&mesh,E,nu);
   ScalarView viewMu("mu [Pa]", new WinGeom(550, 460, 530, 350));
   viewMu.show_mesh(false);
   viewMu.show(&exact_mu);

   // Initialize boundary conditions.
   DefaultEssentialBCConst<double> disp_bot_top_x(Hermes::vector<std::string>("Bottom","Top"), 0.0);
   DefaultEssentialBCConst<double> disp_bot_y("Bottom", 0.0);
   DefaultEssentialBCConst<double> disp_top_y("Top", 0.1);
   EssentialBCs<double> bcs_x(&disp_bot_top_x);
   EssentialBCs<double> bcs_y(Hermes::vector<EssentialBoundaryCondition<double> *>(&disp_bot_y, &disp_top_y));

   // Create x- and y- displacement space using the default H1 shapeset.
   H1Space<double> u1_space(&mesh, &bcs_x, P_INIT);
   H1Space<double> u2_space(&mesh, &bcs_y, P_INIT);
   int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u1_space, &u2_space));
   Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

   // Initialize the weak formulation.
   CustomWeakFormLinearElasticity wf(E, nu, &lambdaFun, &muFun, rho*g1, "Top", f0, f1);

   // Initialize the FE problem.
   DiscreteProblem<double> dp(&wf, Hermes::vector<const Space<double> *>(&u1_space, &u2_space));

   // Initialize Newton solver.
   NewtonSolver<double> newton(&dp);
   newton.set_verbose_output(true);

   // Perform Newton's iteration.
   try
   {
     // NULL = start from zero initial vector, 1e-7 = tolerance.
//     newton.solve(NULL, 1e-7);
	   newton.solve();
   }
//   catch(Hermes::Exceptions::Exception e)
   catch(std::exception& e)
   {
     std::cout << e.what();
	 Hermes::Mixins::Loggable::Static::info("Newton's iteration failed.");
   }

   // Translate the resulting coefficient vector into the Solution sln.
   Solution<double> u1_sln, u2_sln;
   Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<const Space<double> *>(&u1_space, &u2_space),
       Hermes::vector<Solution<double> *>(&u1_sln, &u2_sln));

   // Visualize the solution.
   ScalarView view("Von Mises stress [Pa]", new WinGeom(590, 0, 700, 400));
   // First Lame constant.
   double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
   // Second Lame constant.
   double mu = E / (2*(1 + nu));
   VonMisesFilter stress(Hermes::vector<MeshFunction<double> *>(&u1_sln, &u2_sln), lambda, mu);
   view.show_mesh(false);
   view.show(&stress, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u1_sln, &u2_sln, 1.0);

   ScalarView viewS11("S11 [Pa]", new WinGeom(0, 260, 530, 350));
   CustomFilterS11 S11(Hermes::vector<Solution<double> *>(&u1_sln, &u2_sln), mu, lambda);
   viewS11.show_mesh(false);
   viewS11.show(&S11, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u1_sln, &u2_sln, 1.0);

   ScalarView viewS12("S12 [Pa]", new WinGeom(540, 260, 530, 350));
   CustomFilterS12 S12(Hermes::vector<Solution<double> *>(&u1_sln, &u2_sln), mu);
   viewS12.show_mesh(false);
   viewS12.show(&S12, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u1_sln, &u2_sln, 1.0);

   ScalarView viewS22("S22 [Pa]", new WinGeom(1080, 260, 530, 350));
   CustomFilterS22 S22(Hermes::vector<Solution<double> *>(&u1_sln, &u2_sln), mu, lambda);
   viewS22.show_mesh(false);
   viewS22.show(&S22, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u1_sln, &u2_sln, 1.0);

   // Wait for the view to be closed.
   View::wait();

   return 0;
}
