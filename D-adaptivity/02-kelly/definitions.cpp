#include "definitions.h"

CustomWeakFormPoisson::CustomWeakFormPoisson(const std::string& mat_motor, double eps_motor,
  const std::string& mat_air, double eps_air, Mesh* mesh)
  : WeakForm<double>(1), mat_motor(mat_motor), eps_motor(eps_motor), mat_air(mat_air), eps_air(eps_air), mesh(mesh)
{
    // Jacobian.
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, mat_motor, new Hermes1DFunction<double>(eps_motor)));
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, mat_air, new Hermes1DFunction<double>(eps_air)));

    // Residual.
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, mat_motor, new Hermes1DFunction<double>(eps_motor)));
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, mat_air, new Hermes1DFunction<double>(eps_air)));
  }

double CustomWeakFormPoisson::get_element_eps(Hermes2D::Geom< double >* e)
{
  std::string marker = mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;

  if (marker == mat_motor)
    return eps_motor;
  else if (marker == mat_air)
    return eps_air;

  error("Unknown element marker.");
  return -1;
}

double ResidualErrorForm::value(int n, double* wt,
  Func< double >* u_ext[], Func< double >* u,
  GeomVol< double >* e, ExtData< double >* ext) const
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  double result = 0.;

  // Calculate the integral of squared residual: (f + eps*laplace u)^2.
  for (int i = 0; i < n; i++)
    result += wt[i] * sqr(0.0 + eps * u->laplace[i]);

  return result * sqr(e->diam) / 24.;
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in hermes2d_common_defs.h"
    "if you want to use second derivatives in weak forms.");
#endif
}

Ord ResidualErrorForm::ord(int n, double* wt,
  Func< Ord >* u_ext[], Func< Ord >* u,
  GeomVol< Ord >* e, ExtData< Ord >* ext) const
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  return sqr(u->laplace[0]);
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in hermes2d_common_defs.h"
    "if you want to use second derivatives in weak forms.");
#endif
}

double EnergyErrorForm::value(int n, double* wt,
  Func< double >* u_ext[], Func< double >* u, Func< double >* v,
  GeomVol< double >* e, ExtData< double >* ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * wf->get_element_eps(e) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);

  return result;
}

Ord EnergyErrorForm::ord(int n, double* wt,
  Func< Ord >* u_ext[], Func< Ord >* u, Func< Ord >* v,
  GeomVol< Ord >* e, ExtData< Ord >* ext) const
{
  return u->dx[0] * v->dx[0] + u->dy[0] * v->dy[0];
}