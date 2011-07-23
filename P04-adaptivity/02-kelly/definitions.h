#include "hermes2d.h"

using namespace Hermes::Hermes2D;

class CustomWeakFormPoisson : public WeakForm<double>
{
  public:
    CustomWeakFormPoisson(const std::string& mat_motor, double eps_motor, 
                          const std::string& mat_air, double eps_air);
                          
    double get_element_eps(Geom<double> *e);
    
  private:
    std::string mat_motor;
    std::string mat_air;
    double eps_motor;
    double eps_air;                          
};

class CustomInterfaceEstimatorScalingFunction : public InterfaceEstimatorScalingFunction
{
  public:
    CustomInterfaceEstimatorScalingFunction(const std::string& mat_motor, double eps_motor,
                                            const std::string& mat_air, double eps_air)
      : use_eps(true), mat_motor(mat_motor), eps_motor(eps_motor), mat_air(mat_air), eps_air(eps_air)
    { };
    
    CustomInterfaceEstimatorScalingFunction() : use_eps(false) { };
    
    virtual double value(double e_diam, const std::string& e_marker) const
    { 
      if (use_eps)
      {
        if (e_marker == mat_motor)
          return Hermes::sqr(eps_motor) * e_diam / 24.;
        else if (e_marker == mat_air)
          return Hermes::sqr(eps_air) * e_diam / 24.;
      }
      else
        return e_diam/24.;
      
      return -1; // Never happens, but leave uncommented to suppress compiler warnings.
    }
    
  private:
    bool use_eps;
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

/* Bilinear form inducing the energy norm */

class EnergyErrorForm : public Adapt<double>::MatrixFormVolError
{
public:
  EnergyErrorForm(CustomWeakFormPoisson *wf) 
    : Adapt<double>::MatrixFormVolError(HERMES_UNSET_NORM), wf(wf)
  { };

  virtual double value(int n, double *wt, 
                       Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                       Geom<double> *e, ExtData<double> *ext) const;
  inline virtual Hermes::Ord ord(int n, double *wt, 
                                 Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
                                 Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;
private:
  CustomWeakFormPoisson *wf;
};
