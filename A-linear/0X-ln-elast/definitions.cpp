#include "definitions.h"

CustomWeakFormLinearElasticity::CustomWeakFormLinearElasticity(double E, double nu, Hermes2DFunction<double>* lambdaF, Hermes2DFunction<double>* muF, double rho_g,
                                 std::string surface_force_bdy, double f0, double f1) : WeakForm<double>(2)
{
  double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
  double mu = E / (2*(1 + nu));

  // SINGLE-COMPONENT FORMS. USEFUL FOR MULTIMESH, DO NOT REMOVE.
  // Jacobian.
  add_matrix_form(new CustomJacobianElast00(0, 0, lambda, mu, lambdaF, muF));
  add_matrix_form(new CustomJacobianElast01(0, 1, lambda, mu, lambdaF, muF));
  add_matrix_form(new CustomJacobianElast10(1, 0, lambda, mu, lambdaF, muF));
  add_matrix_form(new CustomJacobianElast11(1, 1, lambda, mu, lambdaF, muF));

  //Residuals
  add_vector_form(new CustomVectorRes0(0, lambda, mu, lambdaF, muF));
  add_vector_form(new CustomVectorRes1(1, lambda, mu, lambdaF, muF));
}

// Jacobian lin elast 0-0
double CustomWeakFormLinearElasticity::CustomJacobianElast00::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * ((lambda + 2 * mu) * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
    result += wt[i] * ((lambdaF->value(e->x[i],e->y[i]) + 2 * muF->value(e->x[i],e->y[i])) * u->dx[i] * v->dx[i] + muF->value(e->x[i],e->y[i]) * u->dy[i] * v->dy[i]);
  }
  return result;
}

Ord CustomWeakFormLinearElasticity::CustomJacobianElast00::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * ((lambda + 2 * mu) * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
	result += wt[i] * ((lambdaF->value(e->x[i],e->y[i]) + 2 * muF->value(e->x[i],e->y[i])) * u->dx[i] * v->dx[i] + muF->value(e->x[i],e->y[i]) * u->dy[i] * v->dy[i]);
  }
  return result;
}

MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast00::clone() const
{
	return new CustomJacobianElast00(*this);
}

// Jacobian lin elast 0-1
double CustomWeakFormLinearElasticity::CustomJacobianElast01::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
	result += wt[i] * (muF->value(e->x[i],e->y[i]) * u->dx[i] * v->dy[i] + lambdaF->value(e->x[i],e->y[i]) * u->dy[i] * v->dx[i]);
  }
  return result;
}

Ord CustomWeakFormLinearElasticity::CustomJacobianElast01::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    result += wt[i] * (muF->value(e->x[i],e->y[i]) * u->dx[i] * v->dy[i] + lambdaF->value(e->x[i],e->y[i]) * u->dy[i] * v->dx[i]);
  }
  return result;
}

MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast01::clone() const
{
	return new CustomJacobianElast01(*this);
}

// Jacobian lin elast 1-0
double CustomWeakFormLinearElasticity::CustomJacobianElast10::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * (mu * u->dy[i] * v->dx[i] + lambda * u->dx[i] * v->dy[i]);
	result += wt[i] * (muF->value(e->x[i],e->y[i]) * u->dy[i] * v->dx[i] + lambdaF->value(e->x[i],e->y[i]) * u->dx[i] * v->dy[i]);
  }
  return result;
}

Ord CustomWeakFormLinearElasticity::CustomJacobianElast10::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * (mu * u->dy[i] * v->dx[i] + lambda * u->dx[i] * v->dy[i]);
    result += wt[i] * (muF->value(e->x[i],e->y[i]) * u->dy[i] * v->dx[i] + lambdaF->value(e->x[i],e->y[i]) * u->dx[i] * v->dy[i]);
  }
  return result;
}

MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast10::clone() const
{
	return new CustomJacobianElast10(*this);
}

// Jacobian lin elast 1-1
double CustomWeakFormLinearElasticity::CustomJacobianElast11::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * ((lambda + 2 * mu) * u->dy[i] * v->dy[i] + mu * u->dx[i] * v->dx[i]);
	result += wt[i] * ((lambdaF->value(e->x[i],e->y[i]) + 2 * muF->value(e->x[i],e->y[i])) * u->dy[i] * v->dy[i] + muF->value(e->x[i],e->y[i]) * u->dx[i] * v->dx[i]);
  }
  return result;
}

Ord CustomWeakFormLinearElasticity::CustomJacobianElast11::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * ((lambda + 2 * mu) * u->dy[i] * v->dy[i] + mu * u->dx[i] * v->dx[i]);
    result += wt[i] * ((lambdaF->value(e->x[i],e->y[i]) + 2 * muF->value(e->x[i],e->y[i])) * u->dy[i] * v->dy[i] + muF->value(e->x[i],e->y[i]) * u->dx[i] * v->dx[i]);
  }
  return result;
}

MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast11::clone() const
{
	return new CustomJacobianElast11(*this);
}

//residuum lin elast 0
double CustomWeakFormLinearElasticity::CustomVectorRes0::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * ((lambda + 2 * mu) * u_ext[0]->dx[i] * v->dx[i] + lambda * u_ext[1]->dy[i] * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i]) * v->dy[i]);
    result += wt[i] * ((lambdaF->value(e->x[i],e->y[i]) + 2 * muF->value(e->x[i],e->y[i])) * u_ext[0]->dx[i] * v->dx[i] + lambdaF->value(e->x[i],e->y[i]) * u_ext[1]->dy[i] * v->dx[i] + muF->value(e->x[i],e->y[i]) * (u_ext[0]->dy[i] + u_ext[1]->dx[i]) * v->dy[i]);
  }
  return result;
}

Ord CustomWeakFormLinearElasticity::CustomVectorRes0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * ((lambda + 2 * mu) * u_ext[0]->dx[i] * v->dx[i] + lambda * u_ext[1]->dy[i] * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i]) * v->dy[i]);
	result += wt[i] * ((lambdaF->value(e->x[i],e->y[i]) + 2 * muF->value(e->x[i],e->y[i])) * u_ext[0]->dx[i] * v->dx[i] + lambdaF->value(e->x[i],e->y[i]) * u_ext[1]->dy[i] * v->dx[i] + muF->value(e->x[i],e->y[i]) * (u_ext[0]->dy[i] + u_ext[1]->dx[i]) * v->dy[i]);
  }
  return result;
}

VectorFormVol<double>* CustomWeakFormLinearElasticity::CustomVectorRes0::clone() const
{
  return new CustomVectorRes0(*this);
}

//residuum lin elast 1
double CustomWeakFormLinearElasticity::CustomVectorRes1::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * ((lambda + 2 * mu) * u_ext[1]->dy[i] * v->dy[i] + lambda * u_ext[0]->dx[i] * v->dy[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i]) * v->dx[i]);
	result += wt[i] * ((lambdaF->value(e->x[i],e->y[i]) + 2 * muF->value(e->x[i],e->y[i])) * u_ext[1]->dy[i] * v->dy[i] + lambdaF->value(e->x[i],e->y[i]) * u_ext[0]->dx[i] * v->dy[i] + muF->value(e->x[i],e->y[i]) * (u_ext[0]->dy[i] + u_ext[1]->dx[i]) * v->dx[i]);
  }
  return result;
}

Ord CustomWeakFormLinearElasticity::CustomVectorRes1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * ((lambda + 2 * mu) * u_ext[0]->dx[i] * v->dx[i] + lambda * u_ext[1]->dy[i] * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i]) * v->dy[i]);
    result += wt[i] * ((lambdaF->value(e->x[i],e->y[i]) + 2 * muF->value(e->x[i],e->y[i])) * u_ext[1]->dy[i] * v->dy[i] + lambdaF->value(e->x[i],e->y[i]) * u_ext[0]->dx[i] * v->dy[i] + muF->value(e->x[i],e->y[i]) * (u_ext[0]->dy[i] + u_ext[1]->dx[i]) * v->dx[i]);
  }
  return result;
}

VectorFormVol<double>* CustomWeakFormLinearElasticity::CustomVectorRes1::clone() const
{
  return new CustomVectorRes1(*this);
}


// Custom function lambda(E,nu)
double CustomExactFunctionLambda::val(double x, double y)
{
  double scaleNu = 10;
  double scaleE = 10;
  double Eln = E * (1-1.0/scaleE) * y + E/scaleE;
  double nuln = nu * (1.0-1.0/scaleNu) * y + nu/scaleNu;
  return (Eln * nuln) / ((1 + nuln) * (1 - 2*nuln));
}
double CustomExactFunctionLambda::dx(double x, double y)
{
  return 0.0;
}
double CustomExactFunctionLambda::ddxx(double x, double y)
{
  return 0.0;
}

// Custom function mu(E,nu)
double CustomExactFunctionMu::val(double x, double y)
{
  double scaleNu = 10;
  double scaleE = 10;
  double Eln = E * (1-1.0/scaleE) * y + E/scaleE;
  double nuln = nu * (1.0-1.0/scaleNu) * y + nu/scaleNu;
  return Eln / (2*(1 + nuln));
}
double CustomExactFunctionMu::dx(double x, double y)
{
  return 0.0;
}
double CustomExactFunctionMu::ddxx(double x, double y)
{
  return 0.0;
}

// Lambda Hermes2DFunction
CustomLambda::CustomLambda(double E, double nu)
  : Hermes2DFunction<double>(), E(E), nu(nu)
{
  cefLam = new CustomExactFunctionLambda(E, nu);
}
double CustomLambda::value(double x, double y) const
{
  return cefLam->val(x,y);
}
Ord CustomLambda::value(Ord x, Ord y) const
{
  return Ord(10);
}
CustomLambda::~CustomLambda()
{
  delete cefLam;
}

// Mu Hermes2DFunction
CustomMu::CustomMu(double E, double nu)
  : Hermes2DFunction<double>(), E(E), nu(nu)
{
  cefMu = new CustomExactFunctionMu(E, nu);
}
double CustomMu::value(double x, double y) const
{
  return cefMu->val(x,y);
}
Ord CustomMu::value(Ord x, Ord y) const
{
  return Ord(10);
}
CustomMu::~CustomMu()
{
  delete cefMu;
}

// Exact solution lambda
ExactSolutionLambda::ExactSolutionLambda(MeshSharedPtr mesh, double E, double nu)
     : ExactSolutionScalar<double>(mesh), E(E), nu(nu)
{
  cefLam = new CustomExactFunctionLambda(E, nu);
}
double ExactSolutionLambda::value(double x, double y) const
{
  return cefLam->val(x,y);
}
void ExactSolutionLambda::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = cefLam->dx(x,y);
  dy = -cefLam->dx(y,y);
}
Ord ExactSolutionLambda::ord(double x, double y) const
{
  return Ord(10);
}
ExactSolutionLambda::~ExactSolutionLambda()
{
  delete cefLam;
}
MeshFunction<double>* ExactSolutionLambda::clone() const
{
  return new ExactSolutionLambda(this->mesh, this->E, this->nu);
}

// Exact solution mu
ExactSolutionMu::ExactSolutionMu(MeshSharedPtr mesh, double E, double nu)
     : ExactSolutionScalar<double>(mesh), E(E), nu(nu)
{
  cefMu = new CustomExactFunctionMu(E, nu);
}
double ExactSolutionMu::value(double x, double y) const
{
  return cefMu->val(x,y);
}
void ExactSolutionMu::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = cefMu->dx(x,y);
  dy = -cefMu->dx(y,y);
}
Ord ExactSolutionMu::ord(double x, double y) const
{
  return Ord(10);
}
ExactSolutionMu::~ExactSolutionMu()
{
  delete cefMu;
}
MeshFunction<double>* ExactSolutionMu::clone() const
{
  return new ExactSolutionMu(this->mesh, this->E, this->nu);
}


// Custom filter S11
void CustomFilterS11::filter_fn(int n, double* x, double* y, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                      double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = 2.0 * this->mu->value(x[i], y[i]) * dx.at(0)[i] + this->lambda->value(x[i], y[i]) * (dx.at(0)[i] + dy.at(1)[i]);
	outdx[i] = 0.0;
	outdy[i] = 0.0;
  }
}

// Custom filter S12
void CustomFilterS12::filter_fn(int n, double* x, double* y, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                      double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
	out[i] = mu * (dy.at(0)[i] + dx.at(1)[i]);
	outdx[i] = 0.0;
	outdy[i] = 0.0;
  }
}

// Custom filter S22
void CustomFilterS22::filter_fn(int n, double* x, double* y, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                      double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
	out[i] = 2.0 * mu * dy.at(1)[i] + lambda * (dy.at(1)[i] + dx.at(0)[i]);
	outdx[i] = 0.0;
	outdy[i] = 0.0;
  }
}
