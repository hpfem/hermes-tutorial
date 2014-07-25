#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(std::vector<std::string> newton_boundaries, double heatcap,
    double rho, double tau, double lambda, double alpha, double temp_ext,
    MeshFunctionSharedPtr<double> sln_prev_time, bool JFNK = false);

  ~CustomWeakForm() {};

private:
  class JacobianFormVol : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol(int i, int j, double heatcap, double rho, double lambda, double tau)
      : MatrixFormVol<double>(i, j), heatcap(heatcap), rho(rho), lambda(lambda), tau(tau) { this->setSymFlag(HERMES_SYM); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const;

    double heatcap, rho, lambda, tau;
  };

  class JacobianFormSurf : public MatrixFormSurf<double>
  {
  public:
    JacobianFormSurf(int i, int j, std::vector<std::string> newton_boundaries, double alpha, double lambda)
      : MatrixFormSurf<double>(i, j), alpha(alpha), lambda(lambda) { this->set_areas(newton_boundaries); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomSurf<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomSurf<Ord> *e, Func<Ord> **ext) const;

    MatrixFormSurf<double>* clone() const;

    double alpha, lambda;
  };

  class ResidualFormVol : public VectorFormVol<double>
  {
  public:
    ResidualFormVol(int i, double heatcap, double rho, double lambda, double tau)
      : VectorFormVol<double>(i), heatcap(heatcap), rho(rho), lambda(lambda), tau(tau)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;

    VectorFormVol<double>* clone() const;

  private:
    double heatcap, rho, lambda, tau;
  };

  class ResidualFormSurf : public VectorFormSurf<double>
  {
  public:
    ResidualFormSurf(int i, std::vector<std::string> newton_boundaries, double alpha, double lambda, double temp_ext)
      : VectorFormSurf<double>(i), alpha(alpha), lambda(lambda), temp_ext(temp_ext)  { this->set_areas(newton_boundaries); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      GeomSurf<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomSurf<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone() const;

  private:
    double alpha, lambda, temp_ext;
  };
};
