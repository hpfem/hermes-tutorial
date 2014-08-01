#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Weak forms */

class CustomWeakFormPoissonDirichlet : public WeakForm<double>
{
public:
  CustomWeakFormPoissonDirichlet(std::string mat_al, double lambda_al,
    std::string mat_cu, double lambda_cu,
    double vol_heat_src);
};

/* Custom non-constant Dirichlet condition */

class CustomDirichletCondition : public EssentialBoundaryCondition<double>
{
public:
  CustomDirichletCondition(std::vector<std::string> markers, double A, double B, double C);

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y) const;

protected:
  double A, B, C;
};
