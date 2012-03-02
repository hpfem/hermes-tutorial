#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// Exact solution.
class CustomExactSolution : public ExactSolutionScalar<double>
{
  public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
  
  virtual MeshFunction<double>* clone()
  {
    return new CustomExactSolution(mesh);
  }
};
