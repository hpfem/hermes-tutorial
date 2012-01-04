#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(Hermes::vector<std::string> newton_boundaries, double heatcap, 
                 double rho, double tau, double lambda, double alpha, double temp_ext, 
                 Solution<double>* sln_prev_time, bool JFNK = false);

  ~CustomWeakForm() {};

private:
  class JacobianFormVol : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol(int i, int j, double heatcap, double rho, double lambda, double tau) 
            : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_SYM), heatcap(heatcap), rho(rho), lambda(lambda), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    MatrixFormVol<double>* clone();

    double heatcap, rho, lambda, tau;
  };

  class JacobianFormSurf : public MatrixFormSurf<double>
  {
  public:
    JacobianFormSurf(int i, int j, Hermes::vector<std::string> newton_boundaries, double alpha, double lambda) 
            : MatrixFormSurf<double>(i, j, newton_boundaries), alpha(alpha), lambda(lambda) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    MatrixFormSurf<double>* clone();

    double alpha, lambda;
  };

  class ResidualFormVol : public VectorFormVol<double>
  {
  public:
    ResidualFormVol(int i, double heatcap, double rho, double lambda, double tau) 
            : VectorFormVol<double>(i), heatcap(heatcap), rho(rho), lambda(lambda), tau(tau)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    VectorFormVol<double>* clone();

  private:
    double heatcap, rho, lambda, tau;
  };
  
  class ResidualFormSurf : public VectorFormSurf<double>
  {
  public:
    ResidualFormSurf(int i, Hermes::vector<std::string> newton_boundaries, double alpha, double lambda, double temp_ext) 
            : VectorFormSurf<double>(i, newton_boundaries), alpha(alpha), lambda(lambda), temp_ext(temp_ext)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    VectorFormSurf<double>* clone();

  private:
    double alpha, lambda, temp_ext;
  };
};

