#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(bool is_matfree = false);
};

