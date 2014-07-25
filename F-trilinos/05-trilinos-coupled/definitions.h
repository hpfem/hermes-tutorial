#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Preconditioners;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(double Le, double alpha, double beta, double kappa, double x1, double tau, bool JFNK, int PRECOND, MeshFunctionSharedPtr<double> omega, MeshFunctionSharedPtr<double> omega_dt, MeshFunctionSharedPtr<double> omega_dc, MeshFunctionSharedPtr<double> t_prev_time_1, MeshFunctionSharedPtr<double> c_prev_time_1, MeshFunctionSharedPtr<double> t_prev_time_2, MeshFunctionSharedPtr<double> c_prev_time_2);

  ~CustomWeakForm() {};

private:
  class PreconditionerForm_0 : public MatrixFormVol<double>
  {
  public:
    PreconditionerForm_0(double tau, double Le)
      : MatrixFormVol<double>(0, 0), tau(tau), Le(Le) { this->setSymFlag(HERMES_SYM); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;

    virtual MatrixFormVol<double>* clone() const { return new PreconditionerForm_0(tau, Le); }
    double tau, Le;
  };

  class PreconditionerForm_1 : public MatrixFormVol<double>
  {
  public:
    PreconditionerForm_1(double tau, double Le)
      : MatrixFormVol<double>(1, 1), tau(tau), Le(Le) { this->setSymFlag(HERMES_SYM); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;
    virtual MatrixFormVol<double>* clone() const { return new PreconditionerForm_1(tau, Le); }

    double tau, Le;
  };

  class JacobianFormVol_0_0 : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol_0_0(double tau)
      : MatrixFormVol<double>(0, 0), tau(tau) { this->setSymFlag(HERMES_SYM); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;
    virtual MatrixFormVol<double>* clone() const { return new JacobianFormVol_0_0(tau); }

    double tau;
  };

  class JacobianFormVol_0_1 : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol_0_1(double tau)
      : MatrixFormVol<double>(0, 1), tau(tau) { this->setSymFlag(HERMES_SYM); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;
    virtual MatrixFormVol<double>* clone() const { return new JacobianFormVol_0_1(tau); }

    double tau;
  };

  class JacobianFormVol_1_0 : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol_1_0(double tau)
      : MatrixFormVol<double>(1, 0), tau(tau) { this->setSymFlag(HERMES_SYM); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;
    virtual MatrixFormVol<double>* clone() const { return new JacobianFormVol_1_0(tau); }

    double tau;
  };

  class JacobianFormVol_1_1 : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol_1_1(double tau, double Le)
      : MatrixFormVol<double>(1, 1), tau(tau), Le(Le) { this->setSymFlag(HERMES_SYM); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;
    virtual MatrixFormVol<double>* clone() const { return new JacobianFormVol_1_1(tau, Le); }

    double tau, Le;
  };

  class JacobianFormSurf_0_0 : public MatrixFormSurf<double>
  {
  public:
    JacobianFormSurf_0_0(std::string bnd_marker, double kappa)
      : MatrixFormSurf<double>(0, 0), kappa(kappa) { this->set_area(bnd_marker); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomSurf<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomSurf<Ord> *e, Func<Ord> **ext) const;
    virtual MatrixFormSurf<double>* clone() const { return new JacobianFormSurf_0_0(this->areas[0], kappa); }

    double kappa;
  };

  class ResidualFormVol_0 : public VectorFormVol<double>
  {
  public:
    ResidualFormVol_0(double tau)
      : VectorFormVol<double>(0), tau(tau)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;
    virtual VectorFormVol<double>* clone() const { return new ResidualFormVol_0(tau); }

  private:
    double tau;
  };

  class ResidualFormVol_1 : public VectorFormVol<double>
  {
  public:
    ResidualFormVol_1(double tau, double Le)
      : VectorFormVol<double>(1), tau(tau), Le(Le)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;
    virtual VectorFormVol<double>* clone() const { return new ResidualFormVol_1(tau, Le); }

  private:
    double tau, Le;
  };

  class ResidualFormSurf_0 : public VectorFormSurf<double>
  {
  public:
    ResidualFormSurf_0(std::string bnd_marker, double kappa)
      : VectorFormSurf<double>(0), kappa(kappa)  { this->set_area(bnd_marker); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      GeomSurf<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomSurf<Ord> *e, Func<Ord> **ext) const;
    virtual VectorFormSurf<double>* clone() const { return new ResidualFormSurf_0(this->areas[0], kappa); }

  private:
    double kappa;
  };

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
};

class InitialSolutionTemperature : public ExactSolutionScalar<double>
{
public:
  InitialSolutionTemperature(MeshSharedPtr mesh, double x1) : ExactSolutionScalar<double>(mesh), x1(x1) {};

  virtual double value(double x, double y) const {
    return (x <= x1) ? 1.0 : exp(x1 - x);
  };

  virtual void derivatives(double x, double y, double& dx, double& dy) const {
    dx = (x <= x1) ? 0.0 : -exp(x1 - x);
    dy = 0.0;
  };

  virtual Ord ord(double x, double y) const {
    return Ord(-exp(x1 - x));
  }

  virtual MeshFunction<double>* clone() const
  {
    return new InitialSolutionTemperature(mesh, x1);
  }

  // Value.
  double x1;
};

class InitialSolutionConcentration : public ExactSolutionScalar<double>
{
public:
  InitialSolutionConcentration(MeshSharedPtr mesh, double x1, double Le) : ExactSolutionScalar<double>(mesh), x1(x1), Le(Le) {};

  virtual double value(double x, double y) const {
    return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x));
  };

  virtual void derivatives(double x, double y, double& dx, double& dy) const {
    dx = (x <= x1) ? 0.0 : Le * exp(x1 - x);
    dy = 0.0;
  };

  virtual Ord ord(double x, double y) const {
    return Ord(exp(Le*(x1 - x)));
  }

  virtual MeshFunction<double>* clone() const
  {
    return new InitialSolutionConcentration(mesh, x1, Le);
  }

  // Value.
  double x1, Le;
};

class CustomFilter : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilter(std::vector<MeshFunctionSharedPtr<double> > solutions, double Le, double alpha, double beta, double kappa, double x1, double tau) : Hermes::Hermes2D::DXDYFilter<double>(solutions), Le(Le), alpha(alpha), beta(beta), kappa(kappa), x1(x1), tau(tau)
  {
  }

  virtual MeshFunction<double>* clone() const
  {
    std::vector<MeshFunctionSharedPtr<double> > slns;
    std::vector<int> items;
    for (int i = 0; i < this->num; i++)
    {
      slns.push_back(dynamic_cast<Solution<double>*>(this->sln[i]->clone()));
    }
    CustomFilter* filter = new CustomFilter(slns, Le, alpha, beta, kappa, x1, tau);
    return filter;
  }

private:
  virtual void filter_fn(int n, double* x, double* y, std::vector<const double *> values, std::vector<const double *> dx, std::vector<const double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
  double tau;
};

class CustomFilterDc : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterDc(std::vector<MeshFunctionSharedPtr<double> > solutions, double Le, double alpha, double beta, double kappa, double x1, double tau) : Hermes::Hermes2D::DXDYFilter<double>(solutions), Le(Le), alpha(alpha), beta(beta), kappa(kappa), x1(x1), tau(tau)
  {
  }

  virtual MeshFunction<double>* clone() const
  {
    std::vector<MeshFunctionSharedPtr<double> > slns;
    std::vector<int> items;
    for (int i = 0; i < this->num; i++)
    {
      slns.push_back(dynamic_cast<Solution<double>*>(this->sln[i]->clone()));
    }
    CustomFilterDc* filter = new CustomFilterDc(slns, Le, alpha, beta, kappa, x1, tau);
    return filter;
  }

private:
  virtual void filter_fn(int n, double* x, double* y, std::vector<const double *> values, std::vector<const double *> dx, std::vector<const double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
  double tau;
};

class CustomFilterDt : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterDt(std::vector<MeshFunctionSharedPtr<double> > solutions, double Le, double alpha, double beta, double kappa, double x1, double tau) : Hermes::Hermes2D::DXDYFilter<double>(solutions), Le(Le), alpha(alpha), beta(beta), kappa(kappa), x1(x1), tau(tau)
  {
  }
  virtual MeshFunction<double>* clone() const
  {
    std::vector<MeshFunctionSharedPtr<double> > slns;
    std::vector<int> items;
    for (int i = 0; i < this->num; i++)
    {
      slns.push_back(dynamic_cast<Solution<double>*>(this->sln[i]->clone()));
    }
    CustomFilterDt* filter = new CustomFilterDt(slns, Le, alpha, beta, kappa, x1, tau);
    return filter;
  }

private:
  virtual void filter_fn(int n, double* x, double* y, std::vector<const double *> values, std::vector<const double *> dx, std::vector<const double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
  double tau;
};