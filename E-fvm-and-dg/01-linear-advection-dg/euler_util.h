#ifndef EULER_UTIL_H
#define EULER_UTIL_H

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

class DiscontinuityDetector
{
public:
  /// Constructor.
  DiscontinuityDetector(std::vector<SpaceSharedPtr<double> > spaces, 
                        std::vector<MeshFunctionSharedPtr<double> > solutions);

  /// Destructor.
   ~DiscontinuityDetector();

  /// Return a reference to the inner structures.
  virtual std::set<int>& get_discontinuous_element_ids() = 0;

protected:
  /// Members.
  std::vector<SpaceSharedPtr<double> > spaces;
  std::vector<MeshFunctionSharedPtr<double> > solutions;
  std::set<int> discontinuous_element_ids;
  MeshSharedPtr mesh;
};

class KrivodonovaDiscontinuityDetector : public DiscontinuityDetector
{
public:
  /// Constructor.
  KrivodonovaDiscontinuityDetector(std::vector<SpaceSharedPtr<double> > spaces, 
                        std::vector<MeshFunctionSharedPtr<double> > solutions);

  /// Destructor.
   ~KrivodonovaDiscontinuityDetector();

  /// Return a reference to the inner structures.
  std::set<int>& get_discontinuous_element_ids();
  std::set<int>& get_discontinuous_element_ids(double threshold);

protected:
  /// Calculates relative (w.r.t. the boundary edge_i of the Element e).
  double calculate_relative_flow_direction(Element* e, int edge_i);

  /// Calculates jumps of all solution components across the edge edge_i of the Element e.
  void calculate_jumps(Element* e, int edge_i, double result[1]);

  /// Calculates h.
  double calculate_h(Element* e, int polynomial_order);

  /// Calculates the norm of the solution on the central element.
  void calculate_norms(Element* e, int edge_i, double result[1]);
};

class KuzminDiscontinuityDetector : public DiscontinuityDetector
{
public:
  /// Constructor.
  KuzminDiscontinuityDetector(std::vector<SpaceSharedPtr<double> > spaces, 
                        std::vector<MeshFunctionSharedPtr<double> > solutions, bool limit_all_orders_independently = false);

  /// Destructor.
   ~KuzminDiscontinuityDetector();

  /// Return a reference to the inner structures.
  std::set<int>& get_discontinuous_element_ids();

  /// Return a reference to the inner structures.
  std::set<int>& get_second_order_discontinuous_element_ids();

  /// Returns info about the method.
  bool get_limit_all_orders_independently();
protected:
  /// Center.
  void find_centroid_values(Hermes::Hermes2D::Element* e, double u_c[1]);
  void find_centroid_derivatives(Hermes::Hermes2D::Element* e, double u_dx_c[1], double u_dy_c[1]);
  void find_second_centroid_derivatives(Hermes::Hermes2D::Element* e, double u_dxx_c[1], double u_dxy_c[1], double u_dyy_c[1]);

  /// Vertices.
  void find_vertex_values(Hermes::Hermes2D::Element* e, double vertex_values[1][4]);
  void find_vertex_derivatives(Hermes::Hermes2D::Element* e, double vertex_derivatives[1][4][2]);

  /// Logic - 1st order.
  void find_u_i_min_max_first_order(Hermes::Hermes2D::Element* e, double u_i_min[1][4], double u_i_max[1][4]);
  void find_alpha_i_first_order(double u_i_min[1][4], double u_i_max[1][4], double u_c[1], double u_i[1][4], double alpha_i[1]);
  void find_alpha_i_first_order_real(Hermes::Hermes2D::Element* e, double u_i[1][4], double u_c[1], double u_dx_c[1], double u_dy_c[1], double alpha_i_real[1]);

  /// Logic - 2nd order.
  void find_u_i_min_max_second_order(Hermes::Hermes2D::Element* e, double u_d_i_min[1][4][2], double u_d_i_max[1][4][2]);
  void find_alpha_i_second_order(double u_d_i_min[1][4][2], double u_d_i_max[1][4][2], double u_dx_c[1], double u_dy_c[1], double u_d_i[1][4][2], double alpha_i[1]);
  void find_alpha_i_second_order_real(Hermes::Hermes2D::Element* e, double u_i[1][4][2], double u_dx_c[1], double u_dy_c[1], double u_dxx_c[1], double u_dxy_c[1], double u_dyy_c[1], double alpha_i_real[1]);

private:
  /// For limiting of second order terms.
  std::set<int> second_order_discontinuous_element_ids;
  bool limit_all_orders_independently;
};

class FluxLimiter
{
public:
  /// Enumeration of types.
  /// Used to pick the proper DiscontinuityDetector.
  enum LimitingType
  {
    Krivodonova,
    Kuzmin
  };
  /// Constructor.
  FluxLimiter(LimitingType type, double* solution_vector, std::vector<SpaceSharedPtr<double> > spaces, bool Kuzmin_limit_all_orders_independently = false);
  FluxLimiter(LimitingType type, double* solution_vector, SpaceSharedPtr<double> space, bool Kuzmin_limit_all_orders_independently = false);

  /// Destructor.
   ~FluxLimiter();

  /// Do the limiting.
  /// With the possibility to also limit the spaces from which the spaces in the constructors are refined.
  virtual void limit_according_to_detector(std::vector<SpaceSharedPtr<double> > coarse_spaces_to_limit = std::vector<SpaceSharedPtr<double> >());
  
  /// For Kuzmin's detector.
  virtual void limit_second_orders_according_to_detector(std::vector<SpaceSharedPtr<double> > coarse_spaces_to_limit = std::vector<SpaceSharedPtr<double> >());
  
  void get_limited_solutions(std::vector<MeshFunctionSharedPtr<double> > solutions_to_limit);
  void get_limited_solution(MeshFunctionSharedPtr<double> solution_to_limit);
protected:
  /// Members.
  double* solution_vector;
  std::vector<SpaceSharedPtr<double> > spaces;
  DiscontinuityDetector* detector;
  std::vector<MeshFunctionSharedPtr<double> > limited_solutions;
};

#endif
