#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::WeakFormsElasticity;
using namespace Hermes::Hermes2D::Views;

class CustomWeakFormLinearElasticity : public WeakForm<double>
{
public:
  CustomWeakFormLinearElasticity(double E, double nu, Hermes2DFunction<double>* lambdaF, Hermes2DFunction<double>* muF, double rho_g,
    std::string surface_force_bdy, double f0, double f1);
private:
  class CustomJacobianElast00 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianElast00(int i, int j, double lambda, double mu, Hermes2DFunction<double>* lambdaF, Hermes2DFunction<double>* muF) : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), lambdaF(lambdaF), muF(muF) {};
    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord>* *ext) const;
    virtual MatrixFormVol<double>* clone() const;
  protected:
    double lambda;
    double mu;
    Hermes2DFunction<double>* lambdaF;
    Hermes2DFunction<double>* muF;
  };
  class CustomJacobianElast01 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianElast01(int i, int j, double lambda, double mu, Hermes2DFunction<double>* lambdaF, Hermes2DFunction<double>* muF) : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), lambdaF(lambdaF), muF(muF) {};
    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord>* *ext) const;
    virtual MatrixFormVol<double>* clone() const;
  protected:
    double lambda;
    double mu;
    Hermes2DFunction<double>* lambdaF;
    Hermes2DFunction<double>* muF;
  };
  class CustomJacobianElast10 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianElast10(int i, int j, double lambda, double mu, Hermes2DFunction<double>* lambdaF, Hermes2DFunction<double>* muF) : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), lambdaF(lambdaF), muF(muF) {};
    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord>* *ext) const;
    virtual MatrixFormVol<double>* clone() const;
  protected:
    double lambda;
    double mu;
    Hermes2DFunction<double>* lambdaF;
    Hermes2DFunction<double>* muF;
  };
  class CustomJacobianElast11 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianElast11(int i, int j, double lambda, double mu, Hermes2DFunction<double>* lambdaF, Hermes2DFunction<double>* muF) : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), lambdaF(lambdaF), muF(muF) {};
    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord>* *ext) const;
    virtual MatrixFormVol<double>* clone() const;
  protected:
    double lambda;
    double mu;
    Hermes2DFunction<double>* lambdaF;
    Hermes2DFunction<double>* muF;
  };
  class CustomVectorRes0 : public VectorFormVol<double>
  {
  public:
    CustomVectorRes0(int i, double lambda, double mu, Hermes2DFunction<double>* lambdaF, Hermes2DFunction<double>* muF)
      : VectorFormVol<double>(i), lambda(lambda), mu(mu), lambdaF(lambdaF), muF(muF)
    {
      }
    virtual double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord>* *ext) const;
    virtual VectorFormVol<double>* clone() const;
  private:
    double lambda;
    double mu;
    Hermes2DFunction<double>* lambdaF;
    Hermes2DFunction<double>* muF;
  };
  class CustomVectorRes1 : public VectorFormVol<double>
  {
  public:
    CustomVectorRes1(int i, double lambda, double mu, Hermes2DFunction<double>* lambdaF, Hermes2DFunction<double>* muF)
      : VectorFormVol<double>(i), lambda(lambda), mu(mu), lambdaF(lambdaF), muF(muF)
    {
      }

    virtual double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord>* *ext) const;
    virtual VectorFormVol<double>* clone() const;
  private:
    double lambda;
    double mu;
    Hermes2DFunction<double>* lambdaF;
    Hermes2DFunction<double>* muF;
  };
};

// Custom function lambda(E,nu)
class CustomExactFunctionLambda
{
public:
  CustomExactFunctionLambda(double E, double nu) : E(E), nu(nu) {};
  double val(double x, double y);
  double dx(double x, double y);
  double ddxx(double x, double y);
protected:
  double E, nu;
};

// Custom function mu(E,nu)
class CustomExactFunctionMu
{
public:
  CustomExactFunctionMu(double E, double nu) : E(E), nu(nu) {};
  double val(double x, double y);
  double dx(double x, double y);
  double ddxx(double x, double y);
protected:
  double E, nu;
};

// Lambda Hermes2DFunction
class CustomLambda : public Hermes2DFunction<double>
{
public:
  CustomLambda(double E, double nu);
  virtual double value(double x, double y) const;
  virtual Ord value(Ord x, Ord y) const;
  ~CustomLambda();
  CustomExactFunctionLambda* cefLam;
protected:
  double E, nu;
};

// Mu Hermes2DFunction
class CustomMu : public Hermes2DFunction<double>
{
public:
  CustomMu(double E, double nu);
  virtual double value(double x, double y) const;
  virtual Ord value(Ord x, Ord y) const;
  ~CustomMu();
  CustomExactFunctionMu* cefMu;
protected:
  double E, nu;
};

// Exact solution lambda
class ExactSolutionLambda : public ExactSolutionScalar<double>
{
public:
  ExactSolutionLambda(MeshSharedPtr mesh, double E, double nu);
  virtual double value(double x, double y) const;
  virtual void derivatives(double x, double y, double& dx, double& dy) const;
  virtual Ord ord(double x, double y) const;
  ~ExactSolutionLambda();
  virtual MeshFunction<double>* clone() const;
  CustomExactFunctionLambda* cefLam;
protected:
  double E, nu;
};

// Exact solution mu
class ExactSolutionMu : public ExactSolutionScalar<double>
{
public:
  ExactSolutionMu(MeshSharedPtr mesh, double E, double nu);
  virtual double value(double x, double y) const;
  virtual void derivatives(double x, double y, double& dx, double& dy) const;
  virtual Ord ord(double x, double y) const;
  ~ExactSolutionMu();
  virtual MeshFunction<double>* clone() const;
  CustomExactFunctionMu* cefMu;
protected:
  double E, nu;
};

// Custom filter S11
class CustomFilterS11 : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterS11(std::vector<MeshFunctionSharedPtr<double> > solutions, CustomMu* mu, CustomLambda* lambda) : DXDYFilter(solutions), mu(mu), lambda(lambda)
  {
  }
  virtual MeshFunction<double>* clone() const
  {
    std::vector<MeshFunctionSharedPtr<double> > slns;
    for (int i = 0; i < this->solutions.size(); i++)
      slns.push_back(this->solutions[i]->clone());
    CustomFilterS11* filter = new CustomFilterS11(slns, mu, lambda);
    return filter;
  }
private:
  virtual void filter_fn(int n, double* x, double* y, const std::vector<const double *>& values, const std::vector<const double *>& dx, const std::vector<const double *>& dy, double* rslt, double* rslt_dx, double* rslt_dy);
  CustomMu* mu;
  CustomLambda* lambda;
};

// Custom filter S12
class CustomFilterS12 : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterS12(std::vector<MeshFunctionSharedPtr<double> > solutions, double mu) : DXDYFilter(solutions), mu(mu)
  {
  }
  virtual MeshFunction<double>* clone() const
  {
    std::vector<MeshFunctionSharedPtr<double> > slns;
    for (int i = 0; i < this->solutions.size(); i++)
      slns.push_back(this->solutions[i]->clone());
    CustomFilterS12* filter = new CustomFilterS12(slns, mu);
    return filter;
  }
private:
  virtual void filter_fn(int n, double* x, double* y, const std::vector<const double *>& values, const std::vector<const double *>& dx, const std::vector<const double *>& dy, double* rslt, double* rslt_dx, double* rslt_dy);
  double mu;
};

// Custom filter S22
class CustomFilterS22 : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterS22(std::vector<MeshFunctionSharedPtr<double> > solutions, double mu, double lambda) : DXDYFilter(solutions), mu(mu), lambda(lambda)
  {
  }
  virtual MeshFunction<double>* clone() const
  {
    std::vector<MeshFunctionSharedPtr<double> > slns;
    for (int i = 0; i < this->solutions.size(); i++)
      slns.push_back(this->solutions[i]->clone());
    CustomFilterS22* filter = new CustomFilterS22(slns, mu, lambda);
    return filter;
  }
private:
  virtual void filter_fn(int n, double* x, double* y, const std::vector<const double *>& values, const std::vector<const double *>& dx, const std::vector<const double *>& dy, double* rslt, double* rslt_dx, double* rslt_dy);
  double mu;
  double lambda;
};
