#include "hermes2d.h"
//#include "weakform/weakform.h"
//#include "integrals/h1.h"
//#include "boundaryconditions/essential_bcs.h"
//#include "weakform_library/weakforms_h1.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Initial condition for the Newton's method */

class CustomInitialCondition : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) 
  {
  };

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(double x, double y) const;
  MeshFunction<double>* clone() const;
};

/* Essential boundary conditions */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition<double>
{
public:
  CustomEssentialBCNonConst(std::string marker) 
           : EssentialBoundaryCondition<double>(Hermes::vector<std::string>()) 
  {
    this->markers.push_back(marker);
  }

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const;
};


