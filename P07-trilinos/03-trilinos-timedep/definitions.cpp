#include "definitions.h"

CustomWeakForm::CustomWeakForm(Hermes::vector<std::string> newton_boundaries, double heatcap, double rho, double tau, 
                               double lambda, double alpha, double temp_ext, Solution<double>* sln_prev_time, bool JFNK) : WeakForm<double>(1, JFNK)
{
  // Jacobian forms - volumetric.
  add_matrix_form(new JacobianFormVol(0, 0, heatcap, rho, lambda, tau));

  // Jacobian forms - surface.
  add_matrix_form_surf(new JacobianFormSurf(0, 0, newton_boundaries, alpha, lambda));

  // Residual forms - volumetric.
  ResidualFormVol* res_form = new ResidualFormVol(0, heatcap, rho, lambda, tau);
  res_form->ext.push_back(sln_prev_time);
  add_vector_form(res_form);

  // Residual forms - surface.
  add_vector_form_surf(new ResidualFormSurf(0, newton_boundaries, alpha, lambda, temp_ext));
}


double CustomWeakForm::JacobianFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                              Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (heatcap * rho * u->val[i] * v->val[i] / tau
                     + lambda * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  return result;
}

Ord CustomWeakForm::JacobianFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                         Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  // Returning the sum of the degrees of the basis and test function plus two.
  return Ord(10);
}


double CustomWeakForm::JacobianFormSurf::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                               Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * alpha * lambda * u->val[i] * v->val[i];
  return result;
}

Ord CustomWeakForm::JacobianFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                          Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  // Returning the sum of the degrees of the basis and test function plus two.
  return Ord(10);
}


double CustomWeakForm::ResidualFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                              Geom<double> *e, ExtData<double> *ext) const 
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (heatcap * rho * (u_ext[0]->val[i] - ext->fn[0]->val[i]) * v->val[i] / tau
                       + lambda * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]));
    return result;
}

Ord CustomWeakForm::ResidualFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                         Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  // Returning the sum of the degrees of the test function and solution plus two.
  return Ord(10);
}


double CustomWeakForm::ResidualFormSurf::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                               Geom<double> *e, ExtData<double> *ext) const 
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * alpha * lambda * (u_ext[0]->val[i] - temp_ext) * v->val[i];
  return result;  
}

Ord CustomWeakForm::ResidualFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                          Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  // Returning the sum of the degrees of the test function and solution plus two.
  return Ord(10);
}

