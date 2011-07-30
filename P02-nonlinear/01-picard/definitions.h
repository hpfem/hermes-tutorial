#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Nonlinearity lambda(u) = pow(u, alpha) */

class CustomNonlinearity : public HermesFunction<double>
{
public:
  CustomNonlinearity(double alpha);

  virtual double value(double u) const;

  virtual Ord value(Ord u) const;

  protected:
    double alpha;
};

/* Weak forms */

// NOTE: The linear problem in the Picard's method is 
//       solved using the Newton's method.

class CustomWeakFormPicard : public WeakForm<double>
{
public:
  CustomWeakFormPicard(Solution<double>* prev_iter_sln, HermesFunction<double>* lambda, HermesFunction<double>* f);

private:
  class CustomJacobian : public MatrixFormVol<double>
  {
  public:
    CustomJacobian(int i, int j, HermesFunction<double>* lambda) : MatrixFormVol(i, j), lambda(lambda) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    protected:
      HermesFunction<double>* lambda;
  };

  class CustomResidual : public VectorFormVol<double>
  {
  public:
    CustomResidual(int i, HermesFunction<double>* lambda, HermesFunction<double>* f) 
      : VectorFormVol<double>(i), lambda(lambda), f(f) 
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[],
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    private:
      HermesFunction<double>* lambda;
      HermesFunction<double>* f;
  };
};

/* Essential boundary conditions */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition<double> 
{
public:
  CustomEssentialBCNonConst(std::string marker) 
           : EssentialBoundaryCondition(Hermes::vector<std::string>()) 
  {
    this->markers.push_back(marker);
  };

  virtual EssentialBoundaryCondition<double>::EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const;
};
