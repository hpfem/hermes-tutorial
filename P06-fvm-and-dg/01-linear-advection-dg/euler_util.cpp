#include "euler_util.h"
#include "limits.h"
#include <limits>

DiscontinuityDetector::DiscontinuityDetector(Hermes::vector<Space<double>*> spaces, 
  Hermes::vector<Solution<double>*> solutions) : spaces(spaces), solutions(solutions)
{
};

DiscontinuityDetector::~DiscontinuityDetector()
{};

KrivodonovaDiscontinuityDetector::KrivodonovaDiscontinuityDetector(Hermes::vector<Space<double>*> spaces, 
  Hermes::vector<Solution<double>*> solutions) : DiscontinuityDetector(spaces, solutions)
{
  // A check that all meshes are the same in the spaces.
  unsigned int mesh0_seq = spaces[0]->get_mesh()->get_seq();
  for(unsigned int i = 0; i < spaces.size(); i++)
    if(spaces[i]->get_mesh()->get_seq() != mesh0_seq)
      error("So far DiscontinuityDetector works only for single mesh.");
  mesh = spaces[0]->get_mesh();
};

KrivodonovaDiscontinuityDetector::~KrivodonovaDiscontinuityDetector()
{};

double KrivodonovaDiscontinuityDetector::calculate_h(Element* e, int polynomial_order)
{
  double h = std::sqrt(std::pow(e->vn[(0 + 1) % e->get_num_surf()]->x - e->vn[0]->x, 2) + std::pow(e->vn[(0 + 1) % e->get_num_surf()]->y - e->vn[0]->y, 2));
  for(int edge_i = 0; edge_i < e->get_num_surf(); edge_i++) {
    double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->y - e->vn[edge_i]->y, 2));
    if(edge_length < h)
      h = edge_length;
  }
  return std::pow(h, (0.5 * (H2D_GET_H_ORDER(spaces[0]->get_element_order(e->id)) 
    + 
    H2D_GET_V_ORDER(spaces[0]->get_element_order(e->id)))
    + 1) / 2);
}

std::set<int>& KrivodonovaDiscontinuityDetector::get_discontinuous_element_ids()
{
  return get_discontinuous_element_ids(1.0);
};

std::set<int>& KrivodonovaDiscontinuityDetector::get_discontinuous_element_ids(double threshold)
{
  Element* e;
  for_all_active_elements(e, mesh)
  {
    bool element_inserted = false;
    for(int edge_i = 0; edge_i < e->get_num_surf() && !element_inserted; edge_i++)
      if(calculate_relative_flow_direction(e, edge_i) < 0 && !e->en[edge_i]->bnd)
      {
        double jumps[1];
        calculate_jumps(e, edge_i, jumps);
        double diameter_indicator = calculate_h(e, spaces[0]->get_element_order(e->id));
        double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->y - e->vn[edge_i]->y, 2));
        double norms[1];
        calculate_norms(e, edge_i, norms);

        // Number of component jumps tested.
        unsigned int component_checked_number = 1;
        for(unsigned int component_i = 0; component_i < component_checked_number; component_i++) {
          if(norms[component_i] < 1E-8)
            continue;
          double discontinuity_detector = jumps[component_i] / (diameter_indicator * edge_length * norms[component_i]);
          if(discontinuity_detector > threshold)
          {
            discontinuous_element_ids.insert(e->id);
            element_inserted = true;
            break;
          }
        }
      }
  }
  return discontinuous_element_ids;
};

double KrivodonovaDiscontinuityDetector::calculate_relative_flow_direction(Element* e, int edge_i)
{
  // Set active element to the two solutions (density_vel_x, density_vel_y).
  solutions[1]->set_active_element(e);
  solutions[2]->set_active_element(e);

  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[1]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 5);
  double3* pt = solutions[1]->get_quad_2d()->get_points(eo);
  int np = solutions[1]->get_quad_2d()->get_num_points(eo);

  Geom<double>* geom = init_geom_surf(solutions[1]->get_refmap(), &surf_pos, eo);
  double3* tan = solutions[1]->get_refmap()->get_tangent(surf_pos.surf_num, eo);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * tan[i][2];

  // Calculate.
  Func<double>* density_vel_x = init_fn(solutions[1], eo);
  Func<double>* density_vel_y = init_fn(solutions[2], eo);

  double result = 0.0;
  for(int point_i = 0; point_i < np; point_i++)
    result += jwt[point_i] * density_vel_x->val[point_i] * geom->nx[point_i] + density_vel_y->val[point_i] * geom->ny[point_i];

  geom->free();
  delete geom;
  delete [] jwt;
  density_vel_x->free_fn();
  density_vel_y->free_fn();
  delete density_vel_x;
  delete density_vel_y;

  return result;
};

void KrivodonovaDiscontinuityDetector::calculate_jumps(Element* e, int edge_i, double result[1])
{
  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[0]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 8);
  double3* pt = solutions[0]->get_quad_2d()->get_points(eo);
  int np = solutions[0]->get_quad_2d()->get_num_points(eo);

  // Initialize the NeighborSearch.
  NeighborSearch<double> ns(e, mesh);
  ns.set_active_edge(edge_i);

  // The values to be returned.
  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;
  result[3] = 0.0;

  // Go through all neighbors.
  for(int neighbor_i = 0; neighbor_i < ns.get_num_neighbors(); neighbor_i++) {
    ns.set_active_segment(neighbor_i);

    // Set active element to the solutions.
    solutions[0]->set_active_element(e);
    solutions[1]->set_active_element(e);
    solutions[2]->set_active_element(e);
    solutions[3]->set_active_element(e);

    // Push all the necessary transformations.
    for(unsigned int trf_i = 0; trf_i < ns.get_central_n_trans(neighbor_i); trf_i++) {
      solutions[0]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
      solutions[1]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
      solutions[2]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
      solutions[3]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
    }

    Geom<double>* geom = init_geom_surf(solutions[0]->get_refmap(), &surf_pos, eo);
    double3* tan = solutions[0]->get_refmap()->get_tangent(surf_pos.surf_num, eo);
    double* jwt = new double[np];
    for(int i = 0; i < np; i++)
      jwt[i] = pt[i][2] * tan[i][2];

    // Prepare functions on the central element.
    Func<double>* density = init_fn(solutions[0], eo);
    Func<double>* density_vel_x = init_fn(solutions[1], eo);
    Func<double>* density_vel_y = init_fn(solutions[2], eo);
    Func<double>* energy = init_fn(solutions[3], eo);

    // Set neighbor element to the solutions.
    solutions[0]->set_active_element(ns.get_neighb_el());
    solutions[1]->set_active_element(ns.get_neighb_el());
    solutions[2]->set_active_element(ns.get_neighb_el());
    solutions[3]->set_active_element(ns.get_neighb_el());

    // Push all the necessary transformations.
    for(unsigned int trf_i = 0; trf_i < ns.get_neighbor_n_trans(neighbor_i); trf_i++) {
      solutions[0]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
      solutions[1]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
      solutions[2]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
      solutions[3]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
    }

    // Prepare functions on the neighbor element.
    Func<double>* density_neighbor = init_fn(solutions[0], eo);
    Func<double>* density_vel_x_neighbor = init_fn(solutions[1], eo);
    Func<double>* density_vel_y_neighbor = init_fn(solutions[2], eo);
    Func<double>* energy_neighbor = init_fn(solutions[3], eo);

    DiscontinuousFunc<double> density_discontinuous(density, density_neighbor, true);
    DiscontinuousFunc<double> density_vel_x_discontinuous(density_vel_x, density_vel_x_neighbor, true);
    DiscontinuousFunc<double> density_vel_y_discontinuous(density_vel_y, density_vel_y_neighbor, true);
    DiscontinuousFunc<double> energy_discontinuous(energy, energy_neighbor, true);

    for(int point_i = 0; point_i < np; point_i++) {
      result[0] += jwt[point_i] * std::abs(density_discontinuous.get_val_central(point_i) - density_discontinuous.get_val_neighbor(point_i)); 
      result[1] += jwt[point_i] * std::abs(density_vel_x_discontinuous.get_val_central(point_i) - density_vel_x_discontinuous.get_val_neighbor(point_i));
      result[2] += jwt[point_i] * std::abs(density_vel_y_discontinuous.get_val_central(point_i) - density_vel_y_discontinuous.get_val_neighbor(point_i));
      result[3] += jwt[point_i] * std::abs(energy_discontinuous.get_val_central(point_i) - energy_discontinuous.get_val_neighbor(point_i));
    }

    geom->free();
    delete geom;
    delete [] jwt;
    density->free_fn();
    density_vel_x->free_fn();
    density_vel_y->free_fn();
    energy->free_fn();
    density_neighbor->free_fn();
    density_vel_x_neighbor->free_fn();
    density_vel_y_neighbor->free_fn();
    energy_neighbor->free_fn();

    delete density;
    delete density_vel_x;
    delete density_vel_y;
    delete energy;
    delete density_neighbor;
    delete density_vel_x_neighbor;
    delete density_vel_y_neighbor;
    delete energy_neighbor;
  }

  result[0] = std::abs(result[0]);
  result[1] = std::abs(result[1]);
  result[2] = std::abs(result[2]);
  result[3] = std::abs(result[3]);
};

void KrivodonovaDiscontinuityDetector::calculate_norms(Element* e, int edge_i, double result[1])
{
  // Set active element to the solutions.
  solutions[0]->set_active_element(e);
  solutions[1]->set_active_element(e);
  solutions[2]->set_active_element(e);
  solutions[3]->set_active_element(e);

  // The values to be returned.
  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;
  result[3] = 0.0;

  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[0]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 8);
  double3* pt = solutions[0]->get_quad_2d()->get_points(eo);
  int np = solutions[0]->get_quad_2d()->get_num_points(eo);

  Geom<double>* geom = init_geom_surf(solutions[0]->get_refmap(), &surf_pos, eo);
  double3* tan = solutions[0]->get_refmap()->get_tangent(surf_pos.surf_num, eo);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * tan[i][2];

  // Calculate.
  Func<double>* density = init_fn(solutions[0], eo);
  Func<double>* density_vel_x = init_fn(solutions[1], eo);
  Func<double>* density_vel_y = init_fn(solutions[2], eo);
  Func<double>* energy = init_fn(solutions[3], eo);

  for(int point_i = 0; point_i < np; point_i++) {
    result[0] = std::max(result[0], std::abs(density->val[point_i]));
    result[1] = std::max(result[1], std::abs(density_vel_x->val[point_i]));
    result[2] = std::max(result[2], std::abs(density_vel_y->val[point_i]));
    result[3] = std::max(result[3], std::abs(energy->val[point_i]));
  }

  geom->free();
  delete geom;
  delete [] jwt;

  density->free_fn();
  density_vel_x->free_fn();
  density_vel_y->free_fn();
  energy->free_fn();

  delete density;
  delete density_vel_x;
  delete density_vel_y;
  delete energy;
};

KuzminDiscontinuityDetector::KuzminDiscontinuityDetector(Hermes::vector<Space<double>*> spaces, 
  Hermes::vector<Solution<double>*> solutions, bool limit_all_orders_independently) : DiscontinuityDetector(spaces, solutions), limit_all_orders_independently(limit_all_orders_independently)
{
  // A check that all meshes are the same in the spaces.
  unsigned int mesh0_seq = spaces[0]->get_mesh()->get_seq();
  for(unsigned int i = 0; i < spaces.size(); i++)
    if(spaces[i]->get_mesh()->get_seq() != mesh0_seq)
      error("So far DiscontinuityDetector works only for single mesh.");
  mesh = spaces[0]->get_mesh();
};

KuzminDiscontinuityDetector::~KuzminDiscontinuityDetector()
{};

std::set<int>& KuzminDiscontinuityDetector::get_discontinuous_element_ids()
{
  Element* e;

  for_all_active_elements(e, mesh)
  {
    if(!limit_all_orders_independently)
      if(this->second_order_discontinuous_element_ids.find(e->id) == this->second_order_discontinuous_element_ids.end())
        continue;
    if(e->is_triangle())
      error("So far this limiter is implemented just for quads.");
    double u_c[1], u_dx_c[1], u_dy_c[1];
    find_centroid_values(e, u_c);
    find_centroid_derivatives(e, u_dx_c, u_dy_c);

    // Vertex values.
    double u_i[1][4];
    find_vertex_values(e, u_i);

    // Boundaries for alpha_i calculation.
    double u_i_min_first_order[1][4];
    double u_i_max_first_order[1][4];
    for(int i = 0; i < 1; i++)
      for(int j = 0; j < 4; j++)
      {
        u_i_min_first_order[i][j] = std::numeric_limits<double>::infinity();
        u_i_max_first_order[i][j] = -std::numeric_limits<double>::infinity();
      }
    find_u_i_min_max_first_order(e, u_i_min_first_order, u_i_max_first_order);

    // alpha_i calculation.
    double alpha_i_first_order[1];
    find_alpha_i_first_order(u_i_min_first_order, u_i_max_first_order, u_c, u_i, alpha_i_first_order);

    // measure.
    for(unsigned int i = 0; i < 4; i++)
      if(1.0 > alpha_i_first_order[i])
      {
        // check for sanity.
        if(std::abs(u_c[i]) > 1E-12 && (std::abs(u_dx_c[i]) > 1E-12 || std::abs(u_dy_c[i]) > 1E-12))
          discontinuous_element_ids.insert(e->id);
      }
  }

  return discontinuous_element_ids;
}

std::set<int>& KuzminDiscontinuityDetector::get_second_order_discontinuous_element_ids()
{
  Element* e;

  for_all_active_elements(e, mesh)
  {
    if(e->is_triangle())
      error("So far this limiter is implemented just for quads.");
    double u_dx_c[1], u_dy_c[1], u_dxx_c[1], u_dxy_c[1], u_dyy_c[1];
    find_centroid_derivatives(e, u_dx_c, u_dy_c);
    find_second_centroid_derivatives(e, u_dxx_c, u_dxy_c, u_dyy_c);

    // Vertex values.
    double u_d_i[1][4][2];
    find_vertex_derivatives(e, u_d_i);
    
    // Boundaries for alpha_i calculation.
    double u_d_i_min_second_order[1][4][2];
    double u_d_i_max_second_order[1][4][2];
    for(int i = 0; i < 1; i++)
      for(int j = 0; j < 4; j++)
        for(int k = 0; k < 2; k++)
        {
          u_d_i_min_second_order[i][j][k] = std::numeric_limits<double>::infinity();
          u_d_i_max_second_order[i][j][k] = -std::numeric_limits<double>::infinity();
        }
    find_u_i_min_max_second_order(e, u_d_i_min_second_order, u_d_i_max_second_order);

    // alpha_i calculation.
    double alpha_i_second_order[1];
    find_alpha_i_second_order(u_d_i_min_second_order, u_d_i_max_second_order, u_dx_c, u_dy_c, u_d_i, alpha_i_second_order);

    // measure.
    for(unsigned int i = 0; i < 4; i++)
      if(1.0 > alpha_i_second_order[i])
      {
        // check for sanity.
        if((std::abs(u_dx_c[i]) > 1E-12 || std::abs(u_dy_c[i]) > 1E-12) && (std::abs(u_dxx_c[i]) > 1E-12 || std::abs(u_dxy_c[i]) > 1E-12 || std::abs(u_dyy_c[i]) > 1E-12))
          second_order_discontinuous_element_ids.insert(e->id);
      }
  }

  return second_order_discontinuous_element_ids;
}

bool KuzminDiscontinuityDetector::get_limit_all_orders_independently()
{
  return this->limit_all_orders_independently;
}

void KuzminDiscontinuityDetector::find_centroid_values(Hermes::Hermes2D::Element* e, double u_c[1])
{
  double c_x, c_y;
  double c_ref_x, c_ref_y;
  if(e->get_num_surf() == 3)
  {
      c_x = (0.33333333333333333) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x);
      c_y = (0.33333333333333333) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y);
  }
  else
  {
      c_x = (0.25) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x + e->vn[3]->x);
      c_y = (0.25) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y + e->vn[3]->y);
  }

  for(unsigned int i = 0; i < this->solutions.size(); i++)
  {
    solutions[i]->set_active_element(e);
    solutions[i]->get_refmap()->untransform(e, c_x, c_y, c_ref_x, c_ref_y);
    u_c[i] = solutions[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 0);
  }
}

void KuzminDiscontinuityDetector::find_centroid_derivatives(Hermes::Hermes2D::Element* e, double u_dx_c[1], double u_dy_c[1])
{
  double c_x, c_y;
  double c_ref_x, c_ref_y;
  if(e->get_num_surf() == 3)
  {
      c_x = (0.33333333333333333) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x);
      c_y = (0.33333333333333333) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y);
  }
  else
  {
      c_x = (0.25) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x + e->vn[3]->x);
      c_y = (0.25) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y + e->vn[3]->y);
  }

  for(unsigned int i = 0; i < this->solutions.size(); i++)
  {
    solutions[i]->set_active_element(e);
    solutions[i]->get_refmap()->untransform(e, c_x, c_y, c_ref_x, c_ref_y);
    u_dx_c[i] = solutions[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 1);
    u_dy_c[i] = solutions[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 2);
  }
}

void KuzminDiscontinuityDetector::find_second_centroid_derivatives(Hermes::Hermes2D::Element* e, double u_dxx_c[1], double u_dxy_c[1], double u_dyy_c[1])
{
  double c_x, c_y;
  double c_ref_x, c_ref_y;
  if(e->get_num_surf() == 3)
  {
      c_x = (0.33333333333333333) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x);
      c_y = (0.33333333333333333) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y);
  }
  else
  {
      c_x = (0.25) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x + e->vn[3]->x);
      c_y = (0.25) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y + e->vn[3]->y);
  }

  for(unsigned int i = 0; i < this->solutions.size(); i++)
  {
    solutions[i]->set_active_element(e);
    solutions[i]->get_refmap()->untransform(e, c_x, c_y, c_ref_x, c_ref_y);
    u_dxx_c[i] = solutions[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 3);
    u_dyy_c[i] = solutions[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 4);
    u_dxy_c[i] = solutions[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 5);
  }
}

void KuzminDiscontinuityDetector::find_vertex_values(Hermes::Hermes2D::Element* e, double vertex_values[1][4])
{
  double c_ref_x, c_ref_y;
  for(unsigned int i = 0; i < this->solutions.size(); i++)
  {
    for(unsigned int j = 0; j < e->get_num_surf(); j++)
    {
      solutions[i]->get_refmap()->set_active_element(e);
      solutions[i]->get_refmap()->untransform(e, e->vn[j]->x, e->vn[j]->y, c_ref_x, c_ref_y);
      vertex_values[i][j] = solutions[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 0);
    }
  }
}

void KuzminDiscontinuityDetector::find_vertex_derivatives(Hermes::Hermes2D::Element* e, double vertex_derivatives[1][4][2])
{
  double c_ref_x, c_ref_y;
  for(unsigned int i = 0; i < this->solutions.size(); i++)
  {
    for(unsigned int j = 0; j < e->get_num_surf(); j++)
    {
      solutions[i]->get_refmap()->set_active_element(e);
      solutions[i]->get_refmap()->untransform(e, e->vn[j]->x, e->vn[j]->y, c_ref_x, c_ref_y);
      vertex_derivatives[i][j][0] = solutions[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 1);
      vertex_derivatives[i][j][1] = solutions[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 2);
    }
  }
}

void KuzminDiscontinuityDetector::find_u_i_min_max_first_order(Hermes::Hermes2D::Element* e, double u_i_min[1][4], double u_i_max[1][4])
{
  for(unsigned int j = 0; j < e->get_num_surf(); j++)
  {
    Hermes::Hermes2D::NeighborSearch<double> ns(e, mesh);
    if(e->en[j]->bnd)
      continue;
    ns.set_active_edge(j);

    // First beginning neighbors on every edge.
    double u_c[1];
    ns.set_active_segment(0);
    find_centroid_values(ns.get_neighb_el(), u_c);
    for(unsigned int min_i = 0; min_i < 1; min_i++)
      if(u_i_min[min_i][j] > u_c[min_i])
        u_i_min[min_i][j] = u_c[min_i];
    for(unsigned int max_i = 0; max_i < 1; max_i++)
      if(u_i_max[max_i][j] < u_c[max_i])
        u_i_max[max_i][j] = u_c[max_i];

    // Second end neighbors on every edge.
    ns.set_active_segment(ns.get_num_neighbors() - 1);
    find_centroid_values(ns.get_neighb_el(), u_c);
    for(unsigned int min_i = 0; min_i < 1; min_i++)
      if(u_i_min[min_i][(j + 1) % e->get_num_surf()] > u_c[min_i])
        u_i_min[min_i][(j + 1) % e->get_num_surf()] = u_c[min_i];
    for(unsigned int max_i = 0; max_i < 1; max_i++)
      if(u_i_max[max_i][(j + 1) % e->get_num_surf()] < u_c[max_i])
        u_i_max[max_i][(j + 1) % e->get_num_surf()] = u_c[max_i];

    // Now the hard part, neighbors' neighbors.
    /// \todo This is where it fails for triangles, where it is much more complicated to look for elements sharing a vertex.
    ns.set_active_segment(0);
    Hermes::Hermes2D::NeighborSearch<double> ns_1(ns.get_neighb_el(), mesh);
    if(ns.get_neighb_el()->en[(ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_num_surf()]->bnd)
      continue;
    ns_1.set_active_edge((ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_num_surf());
    ns_1.set_active_segment(0);
    find_centroid_values(ns_1.get_neighb_el(), u_c);
    for(unsigned int min_i = 0; min_i < 1; min_i++)
      if(u_i_min[min_i][j] > u_c[min_i])
        u_i_min[min_i][j] = u_c[min_i];
    for(unsigned int max_i = 0; max_i < 1; max_i++)
      if(u_i_max[max_i][j] < u_c[max_i])
        u_i_max[max_i][j] = u_c[max_i];
  }
}

void KuzminDiscontinuityDetector::find_alpha_i_first_order(double u_i_min[1][4], double u_i_max[1][4], double u_c[1], double u_i[1][4], double alpha_i[1])
{
  for(unsigned int sol_i = 0; sol_i < 1; sol_i++)
  {
    alpha_i[sol_i] = 1;
    for(unsigned int vertex_i = 0; vertex_i < 4; vertex_i++)
    {
      // Sanity checks.
      if(std::abs(u_i[sol_i][vertex_i] - u_c[sol_i]) < 1E-6)
        continue;
      if(std::abs((u_i_min[sol_i][vertex_i] - u_c[sol_i]) / u_c[sol_i]) > 10)
        continue;
      if(std::abs((u_i_max[sol_i][vertex_i] - u_c[sol_i]) / u_c[sol_i]) > 10)
        continue;

      if(u_i[sol_i][vertex_i] < u_c[sol_i])
      {
        if((u_i_min[sol_i][vertex_i] - u_c[sol_i]) / (u_i[sol_i][vertex_i] - u_c[sol_i]) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_i_min[sol_i][vertex_i] - u_c[sol_i]) / (u_i[sol_i][vertex_i] - u_c[sol_i]);
      }
      else
      {
        if((u_i_max[sol_i][vertex_i] - u_c[sol_i]) / (u_i[sol_i][vertex_i] - u_c[sol_i]) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_i_max[sol_i][vertex_i] - u_c[sol_i]) / (u_i[sol_i][vertex_i] - u_c[sol_i]);
      }
    }
  }
}

void KuzminDiscontinuityDetector::find_u_i_min_max_second_order(Hermes::Hermes2D::Element* e, double u_d_i_min[1][4][2], double u_d_i_max[1][4][2])
{
  for(unsigned int j = 0; j < e->get_num_surf(); j++)
  {
    Hermes::Hermes2D::NeighborSearch<double> ns(e, mesh);
    if(e->en[j]->bnd)
      continue;

    ns.set_active_edge(j);

    // First beginning neighbors on every edge.
    double u_dx_c[1], u_dy_c[1];
    ns.set_active_segment(0);
    find_centroid_derivatives(ns.get_neighb_el(), u_dx_c, u_dy_c);
    for(unsigned int min_i = 0; min_i < 1; min_i++)
    {
      if(u_d_i_min[min_i][j][0] > u_dx_c[min_i])
        u_d_i_min[min_i][j][0] = u_dx_c[min_i];
      if(u_d_i_min[min_i][j][1] > u_dy_c[min_i])
        u_d_i_min[min_i][j][1] = u_dy_c[min_i];
    }
    for(unsigned int max_i = 0; max_i < 1; max_i++)
    {
      if(u_d_i_max[max_i][j][0] < u_dx_c[max_i])
        u_d_i_max[max_i][j][0] = u_dx_c[max_i];
      if(u_d_i_max[max_i][j][1] < u_dy_c[max_i])
        u_d_i_max[max_i][j][1] = u_dy_c[max_i];
    }
    // Second end neighbors on every edge.
    ns.set_active_segment(ns.get_num_neighbors() - 1);
    find_centroid_derivatives(ns.get_neighb_el(), u_dx_c, u_dy_c);
    for(unsigned int min_i = 0; min_i < 1; min_i++)
    {
      if(u_d_i_min[min_i][(j + 1) % e->get_num_surf()][0] > u_dx_c[min_i])
        u_d_i_min[min_i][(j + 1) % e->get_num_surf()][0] = u_dx_c[min_i];
      if(u_d_i_min[min_i][(j + 1) % e->get_num_surf()][1] > u_dy_c[min_i])
        u_d_i_min[min_i][(j + 1) % e->get_num_surf()][1] = u_dy_c[min_i];
    }
    for(unsigned int max_i = 0; max_i < 1; max_i++)
    {
      if(u_d_i_max[max_i][(j + 1) % e->get_num_surf()][0] < u_dx_c[max_i])
        u_d_i_max[max_i][(j + 1) % e->get_num_surf()][0] = u_dx_c[max_i];
      if(u_d_i_max[max_i][(j + 1) % e->get_num_surf()][1] < u_dy_c[max_i])
        u_d_i_max[max_i][(j + 1) % e->get_num_surf()][1] = u_dy_c[max_i];
    }

    // Now the hard part, neighbors' neighbors.
    /// \todo This is where it fails for triangles, where it is much more complicated to look for elements sharing a vertex.
    ns.set_active_segment(0);
    Hermes::Hermes2D::NeighborSearch<double> ns_1(ns.get_neighb_el(), mesh);
    if(ns.get_neighb_el()->en[(ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_num_surf()]->bnd)
      continue;
    ns_1.set_active_edge((ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_num_surf());
    ns_1.set_active_segment(0);
    find_centroid_derivatives(ns_1.get_neighb_el(), u_dx_c, u_dy_c);
    for(unsigned int min_i = 0; min_i < 1; min_i++)
    {
      if(u_d_i_min[min_i][j][0] > u_dx_c[min_i])
        u_d_i_min[min_i][j][0] = u_dx_c[min_i];
      if(u_d_i_min[min_i][j][1] > u_dy_c[min_i])
        u_d_i_min[min_i][j][1] = u_dy_c[min_i];
    }
    for(unsigned int max_i = 0; max_i < 1; max_i++)
    {
      if(u_d_i_max[max_i][j][0] < u_dx_c[max_i])
        u_d_i_max[max_i][j][0] = u_dx_c[max_i];
      if(u_d_i_max[max_i][j][1] < u_dy_c[max_i])
        u_d_i_max[max_i][j][1] = u_dy_c[max_i];
    }
  }
}

void KuzminDiscontinuityDetector::find_alpha_i_second_order(double u_d_i_min[1][4][2], double u_d_i_max[1][4][2], double u_dx_c[1], double u_dy_c[1], double u_d_i[1][4][2], double alpha_i[1])
{
  for(unsigned int sol_i = 0; sol_i < 1; sol_i++)
  {
    alpha_i[sol_i] = 1;
    for(unsigned int vertex_i = 0; vertex_i < 4; vertex_i++)
    {
      // Sanity checks.
      if(std::abs(u_dx_c[sol_i]) < 1E-5)
        continue;
      if(std::abs((u_d_i_min[sol_i][vertex_i][0] - u_dx_c[sol_i]) / u_dx_c[sol_i]) > 10)
        continue;
      if(std::abs((u_d_i_max[sol_i][vertex_i][0] - u_dx_c[sol_i]) / u_dx_c[sol_i]) > 10)
        continue;

      if(std::abs(u_dy_c[sol_i]) < 1E-5)
        continue;
      if(std::abs((u_d_i_min[sol_i][vertex_i][1] - u_dy_c[sol_i]) / u_dy_c[sol_i]) > 10)
        continue;
      if(std::abs((u_d_i_max[sol_i][vertex_i][1] - u_dy_c[sol_i]) / u_dy_c[sol_i]) > 10)
        continue;

      // dx.
      if(u_d_i[sol_i][vertex_i][0] < u_dx_c[sol_i])
      {
        if((u_d_i_min[sol_i][vertex_i][0] - u_dx_c[sol_i]) / (u_d_i[sol_i][vertex_i][0] - u_dx_c[sol_i]) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_d_i_min[sol_i][vertex_i][0] - u_dx_c[sol_i]) / (u_d_i[sol_i][vertex_i][0] - u_dx_c[sol_i]);
      }
      else
      {
        if((u_d_i_max[sol_i][vertex_i][0] - u_dx_c[sol_i]) / (u_d_i[sol_i][vertex_i][0] - u_dx_c[sol_i]) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_d_i_max[sol_i][vertex_i][0] - u_dx_c[sol_i]) / (u_d_i[sol_i][vertex_i][0] - u_dx_c[sol_i]);
      }
      // dy.
      if(u_d_i[sol_i][vertex_i][1] < u_dy_c[sol_i])
      {
        if((u_d_i_min[sol_i][vertex_i][1] - u_dy_c[sol_i]) / (u_d_i[sol_i][vertex_i][1] - u_dy_c[sol_i]) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_d_i_min[sol_i][vertex_i][1] - u_dy_c[sol_i]) / (u_d_i[sol_i][vertex_i][1] - u_dy_c[sol_i]);
      }
      else
      {
        if((u_d_i_max[sol_i][vertex_i][1] - u_dy_c[sol_i]) / (u_d_i[sol_i][vertex_i][1] - u_dy_c[sol_i]) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_d_i_max[sol_i][vertex_i][1] - u_dy_c[sol_i]) / (u_d_i[sol_i][vertex_i] [1]- u_dy_c[sol_i]);
      }
    }
  }
}

FluxLimiter::FluxLimiter(FluxLimiter::LimitingType type, double* solution_vector, Hermes::vector<Space<double>*> spaces, bool Kuzmin_limit_all_orders_independently) : solution_vector(solution_vector), spaces(spaces)
{
  for(unsigned int sol_i = 0; sol_i < spaces.size(); sol_i++)
    limited_solutions.push_back(new Hermes::Hermes2D::Solution<double>(spaces[sol_i]->get_mesh()));

  Solution<double>::vector_to_solutions(solution_vector, spaces, limited_solutions);
  switch(type)
  {
    case Krivodonova:
      this->detector = new KrivodonovaDiscontinuityDetector(spaces, limited_solutions);
      break;
    case Kuzmin:
      this->detector = new KuzminDiscontinuityDetector(spaces, limited_solutions, Kuzmin_limit_all_orders_independently);
      break;
  }
};

FluxLimiter::FluxLimiter(FluxLimiter::LimitingType type, double* solution_vector, Space<double>* space, bool Kuzmin_limit_all_orders_independently) : solution_vector(solution_vector)
{
  spaces.push_back(space);

  for(unsigned int sol_i = 0; sol_i < spaces.size(); sol_i++)
    limited_solutions.push_back(new Hermes::Hermes2D::Solution<double>(spaces[sol_i]->get_mesh()));

  Solution<double>::vector_to_solutions(solution_vector, spaces, limited_solutions);
  switch(type)
  {
    case Krivodonova:
      this->detector = new KrivodonovaDiscontinuityDetector(spaces, limited_solutions);
      break;
    case Kuzmin:
      this->detector = new KuzminDiscontinuityDetector(spaces, limited_solutions, Kuzmin_limit_all_orders_independently);
      break;
  }
};

FluxLimiter::~FluxLimiter()
{
  delete detector;
  for(unsigned int sol_i = 0; sol_i < spaces.size(); sol_i++)
    delete limited_solutions[sol_i];
};

void FluxLimiter::get_limited_solutions(Hermes::vector<Solution<double>*> solutions_to_limit)
{
  for(unsigned int i = 0; i < solutions_to_limit.size(); i++)
    solutions_to_limit[i]->copy(this->limited_solutions[i]);
}

void FluxLimiter::get_limited_solution(Solution<double>* solution_to_limit)
{
  if(this->limited_solutions.size() != 1)
    error("Wrong usage of FluxLimiter::get_limited_solution.");
  solution_to_limit->copy(this->limited_solutions[0]);
}

void FluxLimiter::limit_according_to_detector(Hermes::vector<Space<double> *> coarse_spaces_to_limit)
{
  std::set<int> discontinuous_elements = this->detector->get_discontinuous_element_ids();

  // First adjust the solution_vector.
  for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
    for(std::set<int>::iterator it = discontinuous_elements.begin(); it != discontinuous_elements.end(); it++) {
      AsmList<double> al;
      spaces[space_i]->get_element_assembly_list(spaces[space_i]->get_mesh()->get_element(*it), &al);
      for(unsigned int shape_i = 0; shape_i < al.cnt; shape_i++)
        if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.idx[shape_i])) > 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.idx[shape_i])) > 0)
          solution_vector[al.dof[shape_i]] = 0.0;
    }

    // Now adjust the solutions.
    Solution<double>::vector_to_solutions(solution_vector, spaces, limited_solutions);

    if(coarse_spaces_to_limit != Hermes::vector<Space<double>*>()) {
      // Now set the element order to zero.
      Element* e;

      for_all_elements(e, spaces[0]->get_mesh())
        e->visited = false;

      for(std::set<int>::iterator it = discontinuous_elements.begin(); it != discontinuous_elements.end(); it++) {
        AsmList<double> al;
        spaces[0]->get_element_assembly_list(spaces[0]->get_mesh()->get_element(*it), &al);
        for(unsigned int shape_i = 0; shape_i < al.cnt; shape_i++) {
          if(H2D_GET_H_ORDER(spaces[0]->get_shapeset()->get_order(al.idx[shape_i])) > 0 || H2D_GET_V_ORDER(spaces[0]->get_shapeset()->get_order(al.idx[shape_i])) > 0) {
            spaces[0]->get_mesh()->get_element(*it)->visited = true;
            bool all_sons_visited = true;
            for(unsigned int son_i = 0; son_i < 1; son_i++)
              if(!spaces[0]->get_mesh()->get_element(*it)->parent->sons[son_i]->visited)
              {
                all_sons_visited = false;
                break;
              }
              if(all_sons_visited)
                for(unsigned int space_i = 0; space_i < spaces.size(); space_i++) 
                  coarse_spaces_to_limit[space_i]->set_element_order_internal(spaces[space_i]->get_mesh()->get_element(*it)->parent->id, 0);
          }
        }
      }

      Space<double>::assign_dofs(coarse_spaces_to_limit);
    }
};

void FluxLimiter::limit_second_orders_according_to_detector(Hermes::vector<Space<double> *> coarse_spaces_to_limit)
{
  std::set<int> discontinuous_elements;
  if(dynamic_cast<KuzminDiscontinuityDetector*>(this->detector))
    discontinuous_elements = static_cast<KuzminDiscontinuityDetector*>(this->detector)->get_second_order_discontinuous_element_ids();
  else
    error("limit_second_orders_according_to_detector() is to be used only with Kuzmin's vertex based detector.");

  // First adjust the solution_vector.
  for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
    for(std::set<int>::iterator it = discontinuous_elements.begin(); it != discontinuous_elements.end(); it++) {
      AsmList<double> al;
      spaces[space_i]->get_element_assembly_list(spaces[space_i]->get_mesh()->get_element(*it), &al);
      for(unsigned int shape_i = 0; shape_i < al.cnt; shape_i++)
        if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.idx[shape_i])) > 1 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.idx[shape_i])) > 1)
          solution_vector[al.dof[shape_i]] = 0.0;
    }

    // Now adjust the solutions.
    Solution<double>::vector_to_solutions(solution_vector, spaces, limited_solutions);
    if(dynamic_cast<KuzminDiscontinuityDetector*>(this->detector))
    {
      bool Kuzmin_limit_all_orders_independently = dynamic_cast<KuzminDiscontinuityDetector*>(this->detector)->get_limit_all_orders_independently();
      delete detector;
      this->detector = new KuzminDiscontinuityDetector(spaces, limited_solutions, Kuzmin_limit_all_orders_independently);
    }
    else
    {
      delete detector;
      this->detector = new KrivodonovaDiscontinuityDetector(spaces, limited_solutions);
    }

    if(coarse_spaces_to_limit != Hermes::vector<Space<double>*>()) {
      // Now set the element order to zero.
      Element* e;

      for_all_elements(e, spaces[0]->get_mesh())
        e->visited = false;

      for(std::set<int>::iterator it = discontinuous_elements.begin(); it != discontinuous_elements.end(); it++) {
        AsmList<double> al;
        spaces[0]->get_element_assembly_list(spaces[0]->get_mesh()->get_element(*it), &al);
        for(unsigned int shape_i = 0; shape_i < al.cnt; shape_i++) {
          if(H2D_GET_H_ORDER(spaces[0]->get_shapeset()->get_order(al.idx[shape_i])) > 1 || H2D_GET_V_ORDER(spaces[0]->get_shapeset()->get_order(al.idx[shape_i])) > 1) {
            int h_order_to_set = std::min(1, H2D_GET_H_ORDER(spaces[0]->get_shapeset()->get_order(al.idx[shape_i])));
            int v_order_to_set = std::min(1, H2D_GET_V_ORDER(spaces[0]->get_shapeset()->get_order(al.idx[shape_i])));
            spaces[0]->get_mesh()->get_element(*it)->visited = true;
            bool all_sons_visited = true;
            for(unsigned int son_i = 0; son_i < 1; son_i++)
              if(!spaces[0]->get_mesh()->get_element(*it)->parent->sons[son_i]->visited)
              {
                all_sons_visited = false;
                break;
              }
              if(all_sons_visited)
                for(unsigned int space_i = 0; space_i < spaces.size(); space_i++) 
                  coarse_spaces_to_limit[space_i]->set_element_order_internal(spaces[space_i]->get_mesh()->get_element(*it)->parent->id, H2D_MAKE_QUAD_ORDER(h_order_to_set, v_order_to_set));
          }
        }
      }

      Space<double>::assign_dofs(coarse_spaces_to_limit);
    }
};
