#include "definitions.h"

/* Global functions */

double a_11(double x, double y)
{
  if (y > 0) return 1 + x*x + y*y;
  else return 1;
}
double a_22(double x, double y)
{
  if (y > 0) return 1; else return 1 + x*x + y*y;
}
double a_12(double x, double y)
{
  return 1;
}
double a_21(double x, double y)
{
  return 1;
}
double a_1(double x, double y)
{
  return 0.0;
}
double a_2(double x, double y)
{
  return 0.0;
}
double a_0(double x, double y)
{
  return 0.0;
}

/* Custom non-constant Dirichlet condition */

CustomEssentialBCNonConst::CustomEssentialBCNonConst(std::string marker)
: EssentialBoundaryCondition<double>(marker) { }

inline EssentialBoundaryCondition<double>::EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const
{
  return EssentialBoundaryCondition<double>::BC_FUNCTION;
}

double CustomEssentialBCNonConst::value(double x, double y) const
{
  return -Hermes::cos(M_PI*x);
}

/* Weak forms */

CustomWeakFormGeneral::CustomWeakFormGeneral(std::string bdy_vertical) : WeakForm<double>(1)
{
  // Jacobian forms - volumetric.
  add_matrix_form(new MatrixFormVolGeneral(0, 0));

  // Residual forms - volumetric.
  add_vector_form(new VectorFormVolGeneral(0));

  // Residual forms - surface.
  add_vector_form_surf(new VectorFormSurfGeneral(0, bdy_vertical));
}

CustomWeakFormGeneral::MatrixFormVolGeneral::MatrixFormVolGeneral(int i, int j)
: MatrixFormVol<double>(i, j)
{
  this->setSymFlag(HERMES_SYM);
}

double CustomWeakFormGeneral::MatrixFormVolGeneral::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
  Func<double> *v, GeomVol<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    result += (a_11(x, y) * u->dx[i] * v->dx[i] +
      a_12(x, y) * u->dy[i] * v->dx[i] +
      a_21(x, y) * u->dx[i] * v->dy[i] +
      a_22(x, y) * u->dy[i] * v->dy[i] +
      a_1(x, y) * u->dx[i] * v->val[i] +
      a_2(x, y) * u->dy[i] * v->val[i] +
      a_0(x, y) * u->val[i] * v->val[i]) * wt[i];
  }
  return result;
}

Ord CustomWeakFormGeneral::MatrixFormVolGeneral::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
  GeomVol<Ord> *e, Func<Ord> **ext) const
{
  // Returning the sum of the degrees of the basis and test function plus two.
  return u->val[0] * v->val[0] * e->x[0] * e->x[0];
}

MatrixFormVol<double>* CustomWeakFormGeneral::MatrixFormVolGeneral::clone() const
{
  return new MatrixFormVolGeneral(this->i, this->j);
}

CustomWeakFormGeneral::VectorFormVolGeneral::VectorFormVolGeneral(int i) : VectorFormVol<double>(i)
{
}

double CustomWeakFormGeneral::VectorFormVolGeneral::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
  GeomVol<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    result += wt[i] * (a_11(x, y) * u_ext[0]->dx[i] * v->dx[i] +
      a_12(x, y) * u_ext[0]->dy[i] * v->dx[i] +
      a_21(x, y) * u_ext[0]->dx[i] * v->dy[i] +
      a_22(x, y) * u_ext[0]->dy[i] * v->dy[i] +
      a_1(x, y) * u_ext[0]->dx[i] * v->val[i] +
      a_2(x, y) * u_ext[0]->dy[i] * v->val[i] +
      a_0(x, y) * u_ext[0]->val[i] * v->val[i]);
    result -= wt[i] * rhs(e->x[i], e->y[i]) * v->val[i];
  }
  return result;
}

Ord CustomWeakFormGeneral::VectorFormVolGeneral::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
  GeomVol<Ord> *e, Func<Ord> **ext) const
{
  // Returning the sum of the degrees of the test function and solution plus two.
  return u_ext[0]->val[0] * v->val[0] * e->x[0] * e->x[0];
}

VectorFormVol<double>* CustomWeakFormGeneral::VectorFormVolGeneral::clone() const
{
  return new VectorFormVolGeneral(this->i);
}

double CustomWeakFormGeneral::VectorFormVolGeneral::rhs(double x, double y) const
{
  return 1 + x*x + y*y;
}

CustomWeakFormGeneral::VectorFormSurfGeneral::VectorFormSurfGeneral(int i, std::string area)
: VectorFormSurf<double>(i)
{
  this->set_area(area);
}

double CustomWeakFormGeneral::VectorFormSurfGeneral::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
  GeomSurf<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result -= wt[i] * g_N(e->x[i], e->y[i]) * v->val[i];
  return result;
}

Ord CustomWeakFormGeneral::VectorFormSurfGeneral::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
  GeomSurf<Ord> *e, Func<Ord> **ext) const
{
  // Returning the polynomial degree of the test function plus two.
  return v->val[0] * e->x[0] * e->x[0];
}

double CustomWeakFormGeneral::VectorFormSurfGeneral::g_N(double x, double y) const
{
  return 0;
}

VectorFormSurf<double>* CustomWeakFormGeneral::VectorFormSurfGeneral::clone() const
{
  return new VectorFormSurfGeneral(this->i);
}