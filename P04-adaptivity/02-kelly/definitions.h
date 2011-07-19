#include "hermes2d.h"

using namespace Hermes::Hermes2D;

class CustomWeakFormPoisson : public WeakForm<double>
{
  public:
    CustomWeakFormPoisson(const std::string& mat_motor, double eps_motor, 
                          const std::string& mat_air, double eps_air);
};

class CustomInterfaceEstimatorScalingFunction : public InterfaceEstimatorScalingFunction
{
  public:
    CustomInterfaceEstimatorScalingFunction(const std::string& mat_motor, double eps_motor,
                                            const std::string& mat_air, double eps_air)
      : mat_motor(mat_motor), eps_motor(eps_motor), mat_air(mat_air), eps_air(eps_air)
    { };
    
    virtual double value(double e_diam, const std::string& e_marker) const
    {
      /*
      if (e_marker == mat_motor)
        return Hermes::sqr(eps_motor) * e_diam / 24.;
      else if (e_marker == mat_air)
        return Hermes::sqr(eps_air) * e_diam / 24.;*/
      return e_diam/24.;
    }
    
  private:
    std::string mat_motor;
    std::string mat_air;
    double eps_motor;
    double eps_air;
};

/* Linear form for the residual error estimator */

class ResidualErrorForm : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
  public:
    ResidualErrorForm(const std::string& material, double eps) 
    : KellyTypeAdapt<double>::ErrorEstimatorForm(0, material), eps(eps)
    { };
    
    virtual double value(int n, double *wt, 
                         Func<double> *u_ext[], Func<double> *u, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, 
                            Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

  protected:
    double eps;
};

class ResidualErrorFormMotor : public ResidualErrorForm
{
  public:
    ResidualErrorFormMotor(const std::string& mat_motor, double eps_motor) 
      : ResidualErrorForm(mat_motor, eps_motor)
    { };
};

class ResidualErrorFormAir : public ResidualErrorForm
{
  public:
    ResidualErrorFormAir(const std::string& mat_air, double eps_air) 
      : ResidualErrorForm(mat_air, eps_air)
    { };
};
