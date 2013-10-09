#include "definitions.h"

double ExactSolutionCustom::value (double x, double y) const 
{
  return std::pow(x*x + y*y, 0.25);
}

void ExactSolutionCustom::derivatives (double x, double y, double& dx, double& dy) const 
{
  dx = 0.25 * std::pow(x*x + y*y, -0.75) * 2 * x;
  dy = 0.25 * std::pow(x*x + y*y, -0.75) * 2 * y;
}

Ord ExactSolutionCustom::ord(double x, double y) const 
{
  return Hermes::Ord(std::pow(x*x + y*y, 0.25));
}

MeshFunction<double>* ExactSolutionCustom::clone() const
{
  return new ExactSolutionCustom(this->mesh);
}
