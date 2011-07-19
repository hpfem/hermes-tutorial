#include "definitions.h"

CustomWeakFormPoisson::CustomWeakFormPoisson(const std::string& mat_motor, double eps_motor, 
                                             const std::string& mat_air, double eps_air) : WeakForm<double>(1)
{
  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, mat_motor, new HermesFunction<double>(eps_motor)));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, mat_air, new HermesFunction<double>(eps_air)));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, mat_motor, new HermesFunction<double>(eps_motor)));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, mat_air, new HermesFunction<double>(eps_air)));
}

double ResidualErrorForm::value(int n, double* wt, 
                                Func< double >* u_ext[], Func< double >* u, 
                                Geom< double >* e, ExtData< double >* ext) const
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  double result = 0.;

  // Calculate the integral of squared residual: (f + eps*laplace u)^2.
  for (int i = 0; i < n; i++)
    result += wt[i] * Hermes::sqr( 0.0 + eps * u->laplace[i] );

  return result * Hermes::sqr(e->diam) / 24.;
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in hermes2d_common_defs.h"
        "if you want to use second derivatives in weak forms.");
#endif
}

Hermes::Ord ResidualErrorForm::ord(int n, double* wt, 
                                   Func< Hermes::Ord >* u_ext[], Func< Hermes::Ord >* u, 
                                   Geom< Hermes::Ord >* e, ExtData< Hermes::Ord >* ext) const
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  return Hermes::sqr(u->laplace[0]);
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in hermes2d_common_defs.h"
        "if you want to use second derivatives in weak forms.");
#endif
}