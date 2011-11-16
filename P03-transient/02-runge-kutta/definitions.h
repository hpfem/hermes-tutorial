#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Initial condition */

class CustomInitialCondition : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition(Mesh* mesh, double const_value) : ExactSolutionScalar<double>(mesh), 
    const_value(const_value)
  {
  };

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double const_value;
};

/* Weak forms */

class CustomWeakFormHeatRK : public WeakForm<double>
{
public:
  CustomWeakFormHeatRK(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
                       double* current_time_ptr, double temp_init, double t_final);

private:
  // This form is custom since it contains time-dependent exterior temperature.
  class CustomFormResidualSurf : public VectorFormSurf<double>
  {
  private:
      double h;
  public:
    CustomFormResidualSurf(int i, std::string area, double alpha, double rho,
                           double heatcap, double* current_time_ptr, double temp_init, double t_final)
          : VectorFormSurf<double>(i, area), alpha(alpha), rho(rho),
                                     heatcap(heatcap), current_time_ptr(current_time_ptr),
                                     temp_init(temp_init), t_final(t_final) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual VectorFormSurf<double>* clone();

    // Time-dependent exterior temperature.
    template<typename Real>
    Real temp_ext(Real t) const;

    // Members.
    double alpha, rho, heatcap, *current_time_ptr, temp_init, t_final;
  };
};

