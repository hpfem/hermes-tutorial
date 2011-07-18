#include "hermes2d.h"

using namespace Hermes::Hermes2D;

class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(const std::string& mat_motor, double eps_motor, 
                        const std::string& mat_air, double eps_air);
};
