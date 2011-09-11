#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  The purpose of this example is to show how to use Trilinos for nonlinear time-dependent coupled PDE systems.
//  Solved by NOX solver via Newton or JFNK, with or without preconditioning.
//
//  PDE: Flame propagation (same as tutorial example 19-newton-timedep-flame).
//
//  Domain: Same as in tutorial example 19-newton-timedep-flame.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 2;         // Number of initial uniform mesh refinements.
const int P_INIT = 2;               // Initial polynomial degree of all mesh elements.
const bool JFNK = true;             // true = jacobian-free method,
                                    // false = Newton
const int PRECOND = 2;              // Preconditioning by jacobian (1) (less GMRES iterations, more time to create precond)
                                    // or by approximation of jacobian (2) (less time for precond creation, more GMRES iters).
                                    // in case of jfnk,
                                    // default Ifpack proconditioner in case of Newton.
const double TAU = 0.05;            // Time step.
const double T_FINAL = 60.0;        // Time interval length.
const bool TRILINOS_OUTPUT = true;  // Display more details about nonlinear and linear solvers.

// Problem parameters.
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> left_t("Left", 1.0);
  EssentialBCs<double> bcs_t(&left_t);

  DefaultEssentialBCConst<double> left_c("Left", 0.0);
  EssentialBCs<double> bcs_c(&left_c);

  // Create H1 spaces with default shapesets.
  H1Space<double>* t_space = new H1Space<double>(&mesh, &bcs_t, P_INIT);
  H1Space<double>* c_space = new H1Space<double>(&mesh, &bcs_c, P_INIT);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<Space<double>*>(t_space, c_space));
  info("ndof = %d.", ndof);

  // Define initial conditions.
  InitialSolutionTemperature t_prev_time_1(&mesh, x1);
  InitialSolutionConcentration c_prev_time_1(&mesh, x1, Le);
  InitialSolutionTemperature t_prev_time_2(&mesh, x1);
  InitialSolutionConcentration c_prev_time_2(&mesh, x1, Le);
  Solution<double> t_prev_newton;
  Solution<double> c_prev_newton;

  // Filters for the reaction rate omega and its derivatives.
  CustomFilter omega(Hermes::vector<Solution<double>*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);
  CustomFilterDt omega_dt(Hermes::vector<Solution<double>*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);
  CustomFilterDc omega_dc(Hermes::vector<Solution<double>*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);

  // Initialize visualization.
  ScalarView rview("Reaction rate", new WinGeom(0, 0, 800, 230));

  // Initialize weak formulation.
  CustomWeakForm wf(Le, alpha, beta, kappa, x1, TAU, JFNK, PRECOND, &omega, &omega_dt, &omega_dc, &t_prev_time_1, &c_prev_time_1, &t_prev_time_2, &c_prev_time_2);

  // Project the functions "t_prev_time_1" and "c_prev_time_1" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solutions on the FE meshes.");
  double* coeff_vec = new double[ndof];
  OGProjection<double>::project_global(Hermes::vector<Space<double> *>(t_space, c_space), 
                                       Hermes::vector<MeshFunction<double>*>(&t_prev_time_1, &c_prev_time_1),
                                       coeff_vec);

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize finite element problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<Space<double>*>(t_space, c_space));

  // Initialize NOX solver and preconditioner.
  NewtonSolverNOX<double> solver(&dp);
  MlPrecond<double> pc("sa");
  if (PRECOND)
  {
    if (JFNK) solver.set_precond(pc);
    else solver.set_precond("Ifpack");
  }
  if (TRILINOS_OUTPUT)
    solver.set_output_flags(NOX::Utils::Error | NOX::Utils::OuterIteration |
                            NOX::Utils::OuterIterationStatusTest |
                            NOX::Utils::LinearSolverDetails);

  // Time stepping loop:
  double total_time = 0.0;
  cpu_time.tick_reset();
  for (int ts = 1; total_time <= T_FINAL; ts++)
  {
    info("---- Time step %d, t = %g s", ts, total_time + TAU);

    cpu_time.tick(HERMES_SKIP);
    if (solver.solve(coeff_vec))
    {
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), Hermes::vector<Space<double> *>(t_space, c_space), 
                Hermes::vector<Solution<double> *>(&t_prev_newton, &c_prev_newton));

      cpu_time.tick();
      info("Number of nonlin iterations: %d (norm of residual: %g)",
          solver.get_num_iters(), solver.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)",
          solver.get_num_lin_iters(), solver.get_achieved_tol());

      // Time measurement.
      cpu_time.tick(HERMES_SKIP);

			
      // Skip visualization time.
      cpu_time.tick(HERMES_SKIP);

      // Update global time.
      total_time += TAU;

      // Saving solutions for the next time step.
      if(ts > 1)
      {
        t_prev_time_2.copy(&t_prev_time_1);
        c_prev_time_2.copy(&c_prev_time_1);
      }

      t_prev_time_1.copy(&t_prev_newton);
      c_prev_time_1.copy(&c_prev_newton);
      
      // Visualization.
      rview.set_min_max_range(0.0,2.0);
      rview.show(&omega);
      cpu_time.tick(HERMES_SKIP);
    }
    else
      error("NOX failed.");

    info("Total running time for time level %d: %g s.", ts, cpu_time.tick().last());
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

