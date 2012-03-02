#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Weak forms */

class CustomWeakFormPoissonNewton : public WeakForm<double>
{
public:
  CustomWeakFormPoissonNewton(double lambda, double alpha, double T0, std::string bdy_heat_flux);
};
