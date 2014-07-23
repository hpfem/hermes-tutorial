#include "definitions.h"

CustomWeakFormPoissonNewton::CustomWeakFormPoissonNewton(double lambda, double alpha, double T0,
  std::string bdy_heat_flux) : WeakForm<double>(1)
{
    // Jacobian form - volumetric.
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, HERMES_ANY, new Hermes1DFunction<double>(lambda),
      HERMES_SYM, HERMES_AXISYM_Y));

    // Jacobian form - surface.
    add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf<double>(0, 0, bdy_heat_flux, new Hermes2DFunction<double>(alpha),
      HERMES_AXISYM_Y));

    // Residual forms - volumetric.
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, HERMES_ANY, new Hermes1DFunction<double>(lambda),
      HERMES_AXISYM_Y));

    // Residual form - surface.
    add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf<double>(0, bdy_heat_flux, new Hermes2DFunction<double>(alpha),
      HERMES_AXISYM_Y));
    add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf<double>(0, bdy_heat_flux, new Hermes2DFunction<double>(-alpha * T0),
      HERMES_AXISYM_Y));
  };