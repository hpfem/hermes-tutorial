#include "definitions.h"

CustomNonlinearity::CustomNonlinearity(double alpha): Hermes1DFunction<double>()
  {
    this->is_const = false;
    this->alpha = alpha;
  }

double CustomNonlinearity::value(double u) const
  {
    return -1 - std::pow(u, alpha);
  }

Ord CustomNonlinearity::value(Ord u) const
  { 
    return Ord(10);
  }

double CustomNonlinearity::derivative(double u) const
  {
    return -alpha * std::pow(u, alpha - 1.0);
  }

Ord CustomNonlinearity::derivative(Ord u) const
  {
    return Ord(10);
  }


EssentialBCNonConst::EssentialBCNonConst(std::string marker) : EssentialBoundaryCondition<double>(std::vector<std::string>())
  {
    markers.push_back(marker);
  }

EssentialBoundaryCondition<double>::EssentialBCValueType EssentialBCNonConst::get_value_type() const 
  { 
    return EssentialBoundaryCondition<double>::BC_FUNCTION; 
  }

double EssentialBCNonConst::value(double x, double y) const
  {
    return (x + 10) * (y + 10) / 100.;
  }


void CustomInitialCondition::derivatives (double x, double y, double& dx, double& dy) const 
  {
    dx = (y+10)/100.;
    dy = (x+10)/100.;
  };

double CustomInitialCondition::value (double x, double y) const 
  {
    return (x+10)*(y+10)/100.;
  }

Ord CustomInitialCondition::ord(double x, double y) const 
  {
    return Hermes::Ord((x+10)*(y+10)/100.);
  }
  
MeshFunction<double>* CustomInitialCondition::clone() const
{
  return new CustomInitialCondition(this->mesh);
}