#ifndef PCH_HPP
#define PCH_HPP
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iostream>
#include <filesystem>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <list>
#include <CGAL/Arrangement_with_history_2.h>
#include "instance_loader/BasicDataStructures.hpp"
#include "instance_loader/Geometry.hpp"
#include "instance_loader/shapefile_io_operations.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;
typedef Kernel::Segment_2 Segment;
typedef Kernel::Line_2 Line;
typedef Kernel::Vector_2 Vector;
typedef Kernel::Ray_2 Ray;


typedef std::vector<Segment>::iterator Iterator;
typedef CGAL::AABB_segment_primitive_2<Kernel, Iterator> Primitive;
typedef CGAL::AABB_traits_2<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;


template <typename K>
class CustomSegmentTraits : public CGAL::Arr_segment_traits_2<K> {
public:
  typedef CGAL::Arr_segment_traits_2<K> Base;

  class Curve_2 : public Base::Curve_2 {
  public:
    bool frompoly = false;
    int id = -1;
    int source_p = -1;
    int target_p = -1;

    using Base::Curve_2::Curve_2;

    Curve_2() : Base::Curve_2(), frompoly(false), id(-1) {}
    Curve_2(const typename Base::Curve_2& base, bool frompoly, int id)
        : Base::Curve_2(base), frompoly(frompoly), id(id) {}
    Curve_2(const typename Base::Curve_2& base, bool frompoly, int id, int s, int t)
        : Base::Curve_2(base), frompoly(frompoly), id(id), source_p(s), target_p(t) {}

    Curve_2(const typename Base::Point_2& p,
            const typename Base::Point_2& q,
            bool frompoly,
            int id)
        : Base::Curve_2(p, q), frompoly(frompoly), id(id) {}

    Curve_2(const typename Base::Point_2& p,
            const typename Base::Point_2& q,
            bool frompoly, int id, int s, int t)
        : Base::Curve_2(p, q), frompoly(frompoly), id(id), source_p(s), target_p(t) {}

    Curve_2(const Curve_2& other)
        : Base::Curve_2(other), frompoly(other.frompoly), id(other.id),
          source_p(other.source_p), target_p(other.target_p) {}

    Curve_2& operator=(const Curve_2& other) {
      if (this != &other) {
        Base::Curve_2::operator=(other);
        frompoly = other.frompoly;
        id = other.id;
        source_p = other.source_p;
        target_p = other.target_p;
      }
      return *this;
    }
  };
};


typedef CustomSegmentTraits<Kernel> ArrTraits;
typedef CGAL::Arrangement_with_history_2<ArrTraits> Arrangement;
typedef ArrTraits::Curve_2 Curve;




/*
#include <algorithm.hpp>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "instance_loader/BasicDataStructures.hpp"
#include "instance_loader/Geometry.hpp"
#include "instance_loader/shapefile_io_operations.h"


#include <fstream>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;
typedef Kernel::Segment_2 Segment;*/

#endif
