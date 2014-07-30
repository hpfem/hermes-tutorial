#include "definitions.h"

CustomWeakFormHeatRK::CustomWeakFormHeatRK(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
  double* current_time_ptr, double temp_init, double t_final) : WeakForm<double>(1)
{
  // Jacobian volumetric part.
  add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, HERMES_ANY, new Hermes1DFunction<double>(-lambda / (heatcap * rho))));

  // Jacobian surface part.
  add_matrix_form_surf(new DefaultMatrixFormSurf<double>(0, 0, bdy_air, new Hermes2DFunction<double>(-alpha / (heatcap * rho))));

  // Residual - volumetric.
  add_vector_form(new DefaultResidualDiffusion<double>(0, HERMES_ANY, new Hermes1DFunction<double>(-lambda / (heatcap * rho))));

  // Residual - surface.
  add_vector_form_surf(new CustomFormResidualSurf(0, bdy_air, alpha, rho, heatcap,
    current_time_ptr, temp_init, t_final));
}

double CustomWeakFormHeatRK::CustomFormResidualSurf::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomSurf<double> *e,
  Func<double> **ext) const
{
  double T_ext = temp_ext(get_current_stage_time());
  double result = 0;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (T_ext - u_ext[0]->val[i]) * v->val[i];
  }

  return alpha / (rho * heatcap) * result;
}

Ord CustomWeakFormHeatRK::CustomFormResidualSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomSurf<Ord> *e, Func<Ord> **ext) const
{
  Ord T_ext;
  Ord result;

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (T_ext - u_ext[0]->val[i]) * v->val[i];
  }

  return alpha / (rho * heatcap) * result;
}

VectorFormSurf<double>* CustomWeakFormHeatRK::CustomFormResidualSurf::clone() const
{
  return new CustomFormResidualSurf(*this);
}

template<typename Real>
Real CustomWeakFormHeatRK::CustomFormResidualSurf::temp_ext(Real t) const
{
  return temp_init + 10. * Hermes::sin(2 * M_PI*t / t_final);
}