#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(Mesh* mesh, double slope)
             : ExactSolutionScalar<double>(mesh), slope(slope) {};

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double slope;
};

/* Custom function f */

class CustomFunction: public Hermes2DFunction<double>
{
public:
  CustomFunction(double slope)
        : Hermes2DFunction<double>(), slope(slope) {};

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  double slope;
};

