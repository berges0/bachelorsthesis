#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/AABB_traits.h>

#include <vector>
#include <iostream>
#include <optional>
#include "instance_loader/BasicDataStructures.hpp"
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_segment_primitive_2.h>  // Correct for 2D
#include <CGAL/AABB_traits_2.h>             // Correct for 2D

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;
typedef Kernel::Segment_2 Segment;
typedef Kernel::Ray_2 Ray;

typedef std::vector<Segment>::iterator Iterator;
typedef CGAL::AABB_segment_primitive_2<Kernel, Iterator> Primitive;
typedef CGAL::AABB_traits_2<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

class algorithm{
  public:
  algorithm() = default;

  std::vector<Segment> build_segment_vector(const std::vector<SHPLoader::Point2D>& points,
                          const std::vector<std::vector<int>>& polygons);

  std::vector<Segment> ray_shoot_intersection(std::vector<Segment> segments);


};
#endif //ALGORITHM_HPP
