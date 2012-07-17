#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

class ExactSolutionCustom : public ExactSolutionScalar<double>
{
public:
  ExactSolutionCustom(const Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  MeshFunction<double>* clone();
};