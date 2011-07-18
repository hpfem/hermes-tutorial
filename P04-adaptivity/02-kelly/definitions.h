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
/*      if (e_marker == mat_motor)
        return Hermes::sqr(eps_air) * e_diam / 24.;
      else if (e_marker == mat_air)
        return Hermes::sqr(eps_air) * e_diam / 24.;*/
      return e_diam/24.;
    }
    
  private:
    const std::string& mat_motor;
    const std::string& mat_air;
    double eps_motor;
    double eps_air;
};

/* Linear form for the residual error estimator */

class ResidualErrorForm : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
  public:
    ResidualErrorForm(std::string material, double eps) 
    : KellyTypeAdapt<double>::ErrorEstimatorForm(0, material), eps(eps)
    { };
    
    double residual_estimator(int n, double *wt, 
                            Func<double> *u_ext[], Func<double> *u, 
                            Geom<double> *e, ExtData<double> *ext) const;

    Hermes::Ord residual_estimator(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, 
                                  Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[],
                        Func<double> *u, Geom<double> *e,
                        ExtData<double> *ext) const
    {
      return residual_estimator(n, wt, u_ext, u, e, ext);
    }

    virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
                            Func<Hermes::Ord> *u, Geom<Hermes::Ord> *e,
                            ExtData<Hermes::Ord> *ext) const
    {
      return residual_estimator(n, wt, u_ext, u, e, ext);
    }
    
  protected:
    double eps;
};

class ResidualErrorFormMotor : public ResidualErrorForm
{
  public:
    ResidualErrorFormMotor(std::string mat_motor, double eps_motor) 
      : ResidualErrorForm(mat_motor, eps_motor)
    { };
};

class ResidualErrorFormAir : public ResidualErrorForm
{
  public:
    ResidualErrorFormAir(std::string mat_air, double eps_air) 
      : ResidualErrorForm(mat_air, eps_air)
    { };
};