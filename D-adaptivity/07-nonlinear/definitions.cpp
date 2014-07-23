#include "definitions.h"

/* Initial condition */

double InitialSolutionHeatTransfer::value(double x, double y) const
{
  return (x + 10) * (y + 10) / 100. + 2.;
}

void InitialSolutionHeatTransfer::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = (y + 10) / 10.;
  dy = (x + 10) / 10.;
}

Ord InitialSolutionHeatTransfer::ord(double x, double y) const
{
  return Hermes::Ord((x + 10) * (y + 10) / 100. + 2.);
}

MeshFunction<double>* InitialSolutionHeatTransfer::clone() const
{
  return new InitialSolutionHeatTransfer(this->mesh);
}

InitialSolutionHeatTransfer::~InitialSolutionHeatTransfer()
{
}

/* Essential BC */
EssentialBoundaryCondition<double>::EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const
{
  return EssentialBoundaryCondition<double>::BC_FUNCTION;
}

double CustomEssentialBCNonConst::value(double x, double y) const
{
  return (x + 10) * (y + 10) / 100.;
}