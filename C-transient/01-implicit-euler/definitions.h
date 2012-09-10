#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Weak forms */

class CustomWeakFormHeatRK1 : public WeakForm<double>
{
public:
  CustomWeakFormHeatRK1(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
                        double time_step, double* current_time_ptr, double temp_init, double t_final,
                        Solution<double>* prev_time_sln);

private:
  // This form is custom since it contains previous time-level solution.
  class CustomVectorFormVol : public VectorFormVol<double>
  {
  public:
    CustomVectorFormVol(int i, double time_step)
          : VectorFormVol<double>(i), time_step(time_step) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    double time_step;

    VectorFormVol<double>* clone();
  };

  // This form is custom since it contains time-dependent exterior temperature.
  class CustomVectorFormSurf : public VectorFormSurf<double>
  {
  public:
    CustomVectorFormSurf(int i, std::string area, double alpha, double rho, double heatcap,
                         double time_step, double* current_time_ptr, double temp_init, double t_final)
          : VectorFormSurf<double>(i), alpha(alpha), rho(rho), heatcap(heatcap), time_step(time_step), current_time_ptr(current_time_ptr),
                                     temp_init(temp_init), t_final(t_final) 
    {
      this->setArea(area);
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Time-dependent exterior temperature.
    template<typename Real>
    Real temp_ext(Real t) const;

    double alpha, rho, heatcap, time_step, *current_time_ptr, temp_init, t_final;
    
    VectorFormSurf<double>* clone();
  };
};
