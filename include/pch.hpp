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
#include <CGAL/Arrangement_2.h>
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

typedef CGAL::Arr_segment_traits_2<Kernel> ArrTraits;
typedef CGAL::Arrangement_2<ArrTraits> Arrangement;
typedef ArrTraits::X_monotone_curve_2 Curve;




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