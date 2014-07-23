#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Essential boundary conditions */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition<double>
{
public:
  CustomEssentialBCNonConst(std::string marker);

  inline EssentialBCValueType  get_value_type() const;

  virtual double value(double x, double y) const;
};

/* Weak forms */

class CustomWeakFormGeneral : public WeakForm<double>
{
public:
  CustomWeakFormGeneral(std::string bdy_vertical);

private:
  class MatrixFormVolGeneral : public MatrixFormVol<double>
  {
  public:
    MatrixFormVolGeneral(int i, int j);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const;
  };

  class VectorFormVolGeneral : public VectorFormVol<double>
  {
  public:
    VectorFormVolGeneral(int i);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;

    VectorFormVol<double>* clone() const;
  private:
    double rhs(double x, double y) const;
  };

  class VectorFormSurfGeneral : public VectorFormSurf<double>
  {
  public:
    VectorFormSurfGeneral(int i, std::string area = HERMES_ANY);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      GeomSurf<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomSurf<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone() const;
  private:
    double g_N(double x, double y) const;
  };
};
