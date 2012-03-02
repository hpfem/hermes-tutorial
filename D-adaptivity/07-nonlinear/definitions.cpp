#include "definitions.h"

/* Initial condition */

double InitialSolutionHeatTransfer::value (double x, double y) const 
{
  return (x + 10) * (y + 10) / 100. + 2.;
}

void InitialSolutionHeatTransfer::derivatives (double x, double y, double& dx, double& dy) const 
{
  dx = (y + 10) / 10.;
  dy = (x + 10) / 10.;
}

Ord InitialSolutionHeatTransfer::ord(Ord x, Ord y) const 
{
  return (x + 10) * (y + 10) / 100. + 2.;
}

MeshFunction<double>* InitialSolutionHeatTransfer::clone()
{
  Mesh* new_mesh = new Mesh();
  new_mesh->copy(this->mesh);
  return new InitialSolutionHeatTransfer(new_mesh);
}

InitialSolutionHeatTransfer::~InitialSolutionHeatTransfer()
{
  delete mesh;
}

/* Essential BC */
EssentialBoundaryCondition<double>::EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const 
{
  return EssentialBoundaryCondition<double>::BC_FUNCTION;
}

double CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const 
{
  return (x + 10) * (y + 10) / 100.;
}
