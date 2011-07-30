#include "definitions.h"

double CustomInitialCondition::value(double x, double y) const 
{
  return (x+10) * (y+10) / 100. + 2;
}

void CustomInitialCondition::derivatives (double x, double y, double& dx, double& dy) const 
{   
  dx = (y+10) / 100.;
  dy = (x+10) / 100.;
}

Ord CustomInitialCondition::ord(Ord x, Ord y) const 
{
  return x*y;
}

EssentialBoundaryCondition<double>::EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const 
{ 
  return EssentialBoundaryCondition::BC_FUNCTION; 
}

double CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y, 
                                        double t_x, double t_y) const
{
  return (x+10) * (y+10) / 100.;
}
