#include "hermes2d.h"

using namespace Hermes::Hermes2D;

class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string mat_motor, double eps_motor, 
                        std::string mat_air, double eps_air);
};
