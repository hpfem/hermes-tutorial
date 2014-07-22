#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::RefinementSelectors;
using namespace Hermes::Hermes2D::Views;

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(std::string left_bottom_bnd_part, bool DG = true);

private:
  class MatrixFormVol : public Hermes::Hermes2D::MatrixFormVol<double>
  {
  public:
    MatrixFormVol(int i, int j);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, GeomVol<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e, Func<Ord> **ext) const;
    Hermes::Hermes2D::MatrixFormVol<double>* clone() const;
  };

  class VectorFormVol : public Hermes::Hermes2D::VectorFormVol<double>
  {
  public:
    VectorFormVol(int i);

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, GeomVol<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord> **ext) const;
    Hermes::Hermes2D::VectorFormVol<double>* clone() const;

    template<typename Real>
    Real F(Real x, Real y) const;
  };

  class MatrixFormSurface : public Hermes::Hermes2D::MatrixFormSurf<double>
  {
  public:
    MatrixFormSurface(int i, int j);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, GeomSurf<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, GeomSurf<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomSurf<Ord> *e, Func<Ord> **ext) const;
    Hermes::Hermes2D::MatrixFormSurf<double>* clone() const;
  };

  class MatrixFormInterface : public Hermes::Hermes2D::MatrixFormDG<double>
  {
  public:
    MatrixFormInterface(int i, int j);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, DiscontinuousFunc<Scalar> **u_ext, DiscontinuousFunc<Real> *u, DiscontinuousFunc<Real> *v, GeomSurf<Real> *e, DiscontinuousFunc<Scalar> **ext) const;

    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, GeomSurf<double> *e, DiscontinuousFunc<double> **ext) const;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, GeomSurf<Ord> *e, DiscontinuousFunc<Ord> **ext) const;

    MatrixFormDG<double>* clone() const;
  };

  class VectorFormSurface : public Hermes::Hermes2D::VectorFormSurf<double>
  {
  public:
    VectorFormSurface(int i, std::string left_bottom_bnd_part);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomSurf<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomSurf<Ord> *e, Func<Ord> **ext) const;

    template<typename Real>
    Real F(Real x, Real y) const;
    Hermes::Hermes2D::VectorFormSurf<double>* clone() const;
  };

  double calculate_a_dot_v(double x, double y, double vx, double vy) const;

  Ord calculate_a_dot_v(Ord x, Ord y, Ord vx, Ord vy) const;

  double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;

  Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const;
};