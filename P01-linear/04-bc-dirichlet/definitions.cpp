#include "definitions.h"

/* Weak forms */

CustomWeakFormPoissonDirichlet::CustomWeakFormPoissonDirichlet(std::string mat_al, double lambda_al,
                                                               std::string mat_cu, double lambda_cu,
                                                               double vol_heat_src) : WeakForm<double>(1)
{
  // Jacobian forms - volumetric.
  add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, mat_al, new Hermes1DFunction<double>(lambda_al)));
  add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, mat_cu, new Hermes1DFunction<double>(lambda_cu)));

  // Residual forms - volumetric.
  add_vector_form(new DefaultResidualDiffusion<double>(0, mat_al, new Hermes1DFunction<double>(lambda_al)));
  add_vector_form(new DefaultResidualDiffusion<double>(0, mat_cu, new Hermes1DFunction<double>(lambda_cu)));
  add_vector_form(new DefaultVectorFormVol<double>(0, HERMES_ANY, new Hermes2DFunction<double>(-vol_heat_src)));
};

/* Custom non-constant Dirichlet condition */

CustomDirichletCondition::CustomDirichletCondition(Hermes::vector<std::string> markers, 
                                                   double A, double B, double C)
  : EssentialBoundaryCondition(markers), A(A), B(B), C(C) 
{ 
}

EssentialBoundaryCondition<double>::EssentialBCValueType CustomDirichletCondition::get_value_type() const
{ 
  return EssentialBoundaryCondition::BC_FUNCTION; 
}

double CustomDirichletCondition::value(double x, double y, double n_x, double n_y, 
                                       double t_x, double t_y) const 
{
  return A*x + B*y + C;
}
