#include "hermes2d.h"

/* Exact solution */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::RefinementSelectors;
using namespace Hermes::Hermes2D::WeakFormsHcurl;

class CustomExactSolution : public Hermes::Hermes2D::ExactSolutionVector<std::complex<double> >
{
public:
  CustomExactSolution(const Mesh* mesh);
  ~CustomExactSolution();

  virtual Scalar2<std::complex<double> > value(double x, double y) const;

  virtual void derivatives (double x, double y, Scalar2<std::complex<double> >& dx, Scalar2<std::complex<double> >& dy) const;

  virtual Hermes::Ord ord(Hermes::Ord x, Hermes::Ord y) const;
};

/* Weak forms */

class CustomWeakForm : public Hermes::Hermes2D::WeakForm<std::complex<double> >
{
public:
  CustomWeakForm(double mu_r, double kappa);

  class CustomVectorFormSurf : public VectorFormSurf<std::complex<double> >
  {
  public:
    CustomVectorFormSurf();

    virtual std::complex<double> value(int n, double *wt, Func<std::complex<double> > *u_ext[], 
      Func<double> *v, Geom<double> *e, ExtData<std::complex<double> > *ext) const ;

    virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
      Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const ;

    virtual VectorFormSurf<std::complex<double> >* clone();
  };
};
