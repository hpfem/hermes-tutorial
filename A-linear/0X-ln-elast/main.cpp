#include "definitions.h"

// Read original or XML mesh file.
const bool USE_XML_FORMAT = true;
// Initial polynomial degree of all elements.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 5;

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
  MeshSharedPtr mesh(new Mesh), mesh1(new Mesh);
  if (USE_XML_FORMAT == true)
  {
    MeshReaderH2DXML mloader;
    Hermes::Mixins::Loggable::Static::info("Reading mesh in XML format.");
    mloader.load("square.xml", mesh);
  }
  else
  {
    MeshReaderH2D mloader;
    Hermes::Mixins::Loggable::Static::info("Reading mesh in original format.");
    mloader.load("square.mesh", mesh);
  }

  // Perform uniform mesh refinement.
  int refinement_type = 0;
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements(refinement_type);

  // Show mesh.
  MeshView mv("Mesh", new WinGeom(0, 0, 580, 400));
  mv.show(mesh);

  // Exact lambda
  MeshFunctionSharedPtr<double> exact_lambda(new ExactSolutionLambda(mesh, E, nu));
  ScalarView viewLam("lambda [Pa]", new WinGeom(0, 460, 530, 350));
  viewLam.show_mesh(false);
  viewLam.show(exact_lambda);

  // Exact lambda
  MeshFunctionSharedPtr<double> exact_mu(new ExactSolutionMu(mesh, E, nu));
  ScalarView viewMu("mu [Pa]", new WinGeom(550, 460, 530, 350));
  viewMu.show_mesh(false);
  viewMu.show(exact_mu);

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> disp_bot_top_x(std::vector<std::string>({ "Bottom", "Top" }), 0.0);
  DefaultEssentialBCConst<double> disp_bot_y("Bottom", 0.0);
  DefaultEssentialBCConst<double> disp_top_y("Top", 0.1);
  EssentialBCs<double> bcs_x(&disp_bot_top_x);
  EssentialBCs<double> bcs_y({&disp_bot_y, &disp_top_y});

   // Create x- and y- displacement space using the default H1 shapeset.
   SpaceSharedPtr<double> u1_space(new H1Space<double>(mesh, &bcs_x, P_INIT));
   SpaceSharedPtr<double> u2_space(new H1Space<double>(mesh, &bcs_y, P_INIT));
   int ndof = Space<double>::get_num_dofs({ u1_space, u2_space });
   Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

   // Initialize the weak formulation.
   WeakFormSharedPtr<double> wf(new CustomWeakFormLinearElasticity(E, nu, &lambdaFun, &muFun, rho*g1, "Top", f0, f1));

   // Initialize Newton solver.
   NewtonSolver<double> newton(wf, { u1_space, u2_space });
   newton.set_verbose_output(true);

   // Perform Newton's iteration.
   try
   {
	   newton.solve();
   }
   catch(std::exception& e)
   {
     std::cout << e.what();
	 Hermes::Mixins::Loggable::Static::info("Newton's iteration failed.");
   }

   // Translate the resulting coefficient vector into the Solution sln.
   MeshFunctionSharedPtr<double> u1_sln(new Solution<double>), u2_sln(new Solution<double>);

   Solution<double>::vector_to_solutions(newton.get_sln_vector(), { u1_space, u2_space }, { u1_sln, u2_sln });

   // Visualize the solution.
   ScalarView view("Von Mises stress [Pa]", new WinGeom(590, 0, 700, 400));
   // First Lame constant.

   double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
   // Second Lame constant.
   double mu = E / (2*(1 + nu));

   MeshFunctionSharedPtr<double> stress(new VonMisesFilter({ u1_sln, u2_sln }, lambda, mu));
   MeshFunctionSharedPtr<double> S11(new CustomFilterS11({ u1_sln, u2_sln }, &muFun, &lambdaFun));
   MeshFunctionSharedPtr<double> S12(new CustomFilterS12({ u1_sln, u2_sln }, mu));
   MeshFunctionSharedPtr<double> S22(new CustomFilterS22({ u1_sln, u2_sln }, mu, lambda));

   view.show_mesh(false);
   view.set_linearizer_criterion(LinearizerCriterionFixed(3));
   view.show(stress, H2D_FN_VAL_0, u1_sln, u2_sln, 1.0);

   ScalarView viewS11("S11 [Pa]", new WinGeom(0, 260, 530, 350));
   viewS11.show_mesh(false);
   viewS11.set_linearizer_criterion(LinearizerCriterionFixed(2));
   viewS11.show(S11, H2D_FN_VAL_0, u1_sln, u2_sln, 1.0);

   ScalarView viewS12("S12 [Pa]", new WinGeom(540, 260, 530, 350));
   viewS12.show_mesh(false);
   viewS12.set_linearizer_criterion(LinearizerCriterionFixed(2));
   viewS12.show(S12, H2D_FN_VAL_0, u1_sln, u2_sln, 1.0);

   ScalarView viewS22("S22 [Pa]", new WinGeom(1080, 260, 530, 350));
   viewS22.show_mesh(false);
   viewS22.set_linearizer_criterion(LinearizerCriterionFixed(2));
   viewS22.show(S22, H2D_FN_VAL_0, u1_sln, u2_sln, 1.0);

   // Wait for the view to be closed.
   View::wait();

   return 0;
}
