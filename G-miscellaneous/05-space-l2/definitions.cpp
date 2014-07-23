#include "definitions.h"

double CustomExactSolution::value(double x, double y) const
{
  return -std::pow(x, 4) * std::pow(y, 5);
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = 0; // Not needed for L2 projections.
  dy = 0; // Not needed for L2 projections.
}

Ord CustomExactSolution::ord(double x, double y) const
{
  return Hermes::Ord(Hermes::pow(x, 4) * Hermes::pow(y, 5));
}