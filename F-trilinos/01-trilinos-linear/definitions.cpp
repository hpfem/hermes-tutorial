#include "definitions.h"

double CustomExactSolution::value (double x, double y) const 
{
  return x*x + y*y;
}

void CustomExactSolution::derivatives (double x, double y, double& dx, double& dy) const 
{
  dx = 2*x;
  dy = 2*y;
}

Ord CustomExactSolution::ord(double x, double y) const 
{
  return Ord(x*x + y*y);
}


CustomWeakFormPoisson::CustomWeakFormPoisson(bool is_matfree) : WeakForm<double>(1) 
{
  this->is_matfree = is_matfree;

  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(0, HERMES_ANY, new Hermes2DFunction<double>(4.0)));
}


