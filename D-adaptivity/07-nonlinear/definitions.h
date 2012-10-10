#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

/* Initial condition */

class InitialSolutionHeatTransfer : public ExactSolutionScalar<double>
{
public:
  InitialSolutionHeatTransfer(Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};
  ~InitialSolutionHeatTransfer();

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  MeshFunction<double>* clone() const;
};

/* Essential BC */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition<double> {
public:
  CustomEssentialBCNonConst(std::string marker) : EssentialBoundaryCondition<double>(marker) {
  };

  ~CustomEssentialBCNonConst() {};

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;
};
