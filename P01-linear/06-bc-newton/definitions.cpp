#include "definitions.h"

/* Weak forms */

CustomWeakFormPoissonNewton::CustomWeakFormPoissonNewton(std::string mat_al, Hermes1DFunction<double>* lambda_al,
                                                         std::string mat_cu, Hermes1DFunction<double>* lambda_cu,
                                                         Hermes2DFunction<double>* vol_src_term, std::string bdy_heat_flux,
                                                         double alpha, double t_exterior) : WeakForm<double>(1)
{
  // Jacobian forms - volumetric.
  add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, mat_al, lambda_al));
  add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, mat_cu, lambda_cu));

  // Jacobian forms - surface.
  add_matrix_form_surf(new DefaultMatrixFormSurf<double>(0, 0, bdy_heat_flux, new Hermes2DFunction<double>(alpha)));

  // Residual forms - volumetric.
  add_vector_form(new DefaultResidualDiffusion<double>(0, mat_al, lambda_al));
  add_vector_form(new DefaultResidualDiffusion<double>(0, mat_cu, lambda_cu));
  add_vector_form(new DefaultVectorFormVol<double>(0, HERMES_ANY, vol_src_term));

  // Residual forms - surface.
  add_vector_form_surf(new DefaultResidualSurf<double>(0, bdy_heat_flux, new Hermes2DFunction<double>(alpha)));
  add_vector_form_surf(new DefaultVectorFormSurf<double>(0, bdy_heat_flux, new Hermes2DFunction<double>(-alpha * t_exterior)));
};

/* Custom non-constant Dirichlet condition */

CustomDirichletCondition::CustomDirichletCondition(Hermes::vector<std::string> markers, 
                                                   double A, double B, double C)
    : EssentialBoundaryCondition<double>(markers), A(A), B(B), C(C) 
{ 
}

EssentialBoundaryCondition<double>::EssentialBCValueType CustomDirichletCondition::get_value_type() const
{ 
  return EssentialBoundaryCondition<double>::BC_FUNCTION; 
}

double CustomDirichletCondition::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const 
{
  return A*x + B*y + C;
}

