#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Weak forms */

class CustomWeakFormPoissonNeumann : public WeakForm<double>
{
public:
  CustomWeakFormPoissonNeumann(std::string mat_al, Hermes1DFunction<double>* lambda_al,
    std::string mat_cu, Hermes1DFunction<double>* lambda_cu,
    Hermes2DFunction<double>* vol_src_term, std::string bdy_heat_flux,
    Hermes2DFunction<double>* surf_src_term);
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
