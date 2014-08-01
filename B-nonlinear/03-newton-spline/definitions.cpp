#include "definitions.h"

double CustomInitialCondition::value(double x, double y) const
{
  return (x + 10) * (y + 10) / 100. + 2;
}

void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = (y + 10) / 100.;
  dy = (x + 10) / 100.;
}

Ord CustomInitialCondition::ord(double x, double y) const
{
  return Hermes::Ord(x*y);
}

MeshFunction<double>* CustomInitialCondition::clone() const
{
  return new CustomInitialCondition(this->mesh);
}

EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const
{
  return BC_FUNCTION;
}

double CustomEssentialBCNonConst::value(double x, double y) const
{
  return (x + 10) * (y + 10) / 100.;
}