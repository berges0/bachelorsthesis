//
// Created by samuel-berge on 8/10/25.
//

#ifndef PCH_HPP
#define PCH_HPP
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iostream>
#include <filesystem>
#include <string>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <list>
#include <tuple>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <numeric>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_set_2.h>

#include <CGAL/Arr_extended_dcel.h>
#include <typeinfo>
#include <cmath>
#include "instance_loader/BasicDataStructures.hpp"
#include "instance_loader/Geometry.hpp"
#include "instance_loader/shapefile_io_operations.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <stdexcept>
#include <string>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;
typedef Kernel::Segment_2 Segment;
typedef Kernel::Line_2 Line;
typedef Kernel::Vector_2 Vector;
typedef Kernel::Ray_2 Ray;

typedef std::vector<Segment>::iterator Iterator;
typedef CGAL::AABB_segment_primitive_2<Kernel, Iterator> Primitive;
typedef CGAL::AABB_traits_2<Kernel, Primitive> TreeTraits;
typedef CGAL::AABB_tree<TreeTraits> Tree;

struct SkipSegment {
  const Segment* to_skip;

  SkipSegment(const Segment* s) : to_skip(s) {}

  bool operator()(const Primitive::Id& id) const {
    const Segment& s = *id;

    // Skip exact same segment (by pointer)
    if (&s == to_skip) return true;

    // Skip segments that share an endpoint with the segment to skip
    return shares_endpoint(s, *to_skip);
  }

  bool shares_endpoint(const Segment& a, const Segment& b) const {
    return a.source() == b.source() || a.source() == b.target() ||
           a.target() == b.source() || a.target() == b.target();
  }
};

struct Segment_w_info{
    Segment seg;
    bool from_poly;
    int seg_id;
    int poly_id;
    bool shoot_source;
    bool shoot_target;
    Segment_w_info(const Segment& s, bool fp, int seg_id, int poly_id, bool shoot_source,
      bool shoot_target) : seg(s), from_poly(fp), seg_id(seg_id), poly_id(poly_id),  shoot_source(shoot_source),
      shoot_target(shoot_target) {}
};

typedef bool CurveData;              // data attached to original Curve_2
typedef bool XMonotoneCurveData;     // data attached to x-monotone curves

typedef CGAL::Arr_segment_traits_2<Kernel> BaseTraits;
struct Merge {
  XMonotoneCurveData operator()(const XMonotoneCurveData& d1,
                                const XMonotoneCurveData& d2) const {
    std::cout << "Merging data: " << std::endl;
    // for example: choose lower id, or add ids, etc.
    return d1 || d2; // simple merge, just return true if any of the curves has data
  }
};
// Optional: define conversion functor from CurveData → XMonotoneCurveData
struct Convert {
  XMonotoneCurveData operator()(const CurveData& d) const {
    return d;  // trivial if types are same
  }
};

typedef CGAL::Arr_curve_data_traits_2<BaseTraits, XMonotoneCurveData, Merge, CurveData, Convert> Traits;
typedef Traits::Curve_2 Curve;

//define DCEL with face attributes for lookup
// Custom halfedge class with additional data
struct HalfedgeData
{
    bool frompoly;
    bool visited = false;
};

struct VertexData {
    bool is_on_poly = false;
};

struct FaceData {
    int id = -1;
    bool belongs_to_poly = false;
    double area = 0.0; // Optional: area of the face, if needed
    bool visited = false;
    bool added_hole = false; // Flag to indicate if this face has been added as a hole in a polygon
    //Default Constructor
    FaceData(): id(-1), belongs_to_poly(false), area(0) {}
    // Constructor that initializes only the ID
    FaceData(int _id) : id(_id) {}
    FaceData(int _id, bool _belongs_to_poly) : id(_id), belongs_to_poly(_belongs_to_poly) {}
    FaceData(int _id, bool _belongs_to_poly, double _area)
        : id(_id), belongs_to_poly(_belongs_to_poly), area(_area) {}
};

struct Grid {
    double minX, minY, maxX, maxY;
    double stepX=0., stepY=0., lenX, lenY;
    int nr_cells_x=0., nr_cells_y=0., nr_cells=0.;
    double offset_x, offset_y;
    Grid () = default;
    explicit Grid (const std::vector<Segment_w_info> &segments) {
        double min_X = std::numeric_limits<double>::max();
        double max_X = std::numeric_limits<double>::lowest();
        double min_Y = std::numeric_limits<double>::max();
        double max_Y = std::numeric_limits<double>::lowest();

        for (const auto& s : segments) {
            const auto& a = s.seg.source();
            const auto& b = s.seg.target();
            const double ax = CGAL::to_double(a.x()), ay = CGAL::to_double(a.y());
            const double bx = CGAL::to_double(b.x()), by = CGAL::to_double(b.y());
            min_X = std::min(min_X, std::min(ax, bx));
            max_X = std::max(max_X, std::max(ax, bx));
            min_Y = std::min(min_Y, std::min(ay, by));
            max_Y = std::max(max_Y, std::max(ay, by));
        }

        const double len_X = std::max(0.0, max_X - min_X);
        const double len_Y = std::max(0.0, max_Y - min_Y);
        minX = min_X;
        minY = min_Y;
        maxX = max_X;
        maxY = max_Y;
        lenX = len_X;
        lenY = len_Y;
        offset_x = -min_X;
        offset_y = -min_Y;
    }
};

typedef CGAL::Arr_extended_dcel<Traits,
                                VertexData, // data field for vertices
                                HalfedgeData,// data field for halfedge
                                FaceData> //data field for faces
                                                            DCEL;

//typedef CGAL::Arr_face_extended_dcel<Traits, int>           DCEL;
typedef CGAL::Arrangement_with_history_2<Traits, DCEL>       Arrangement;

typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
    boost::no_property,
    boost::property<boost::edge_index_t, std::size_t> > GraphType;

typedef boost::graph_traits<GraphType>::vertex_descriptor VertexDescriptor;
typedef boost::graph_traits<GraphType>::edge_descriptor EdgeDescriptor;
typedef boost::graph_traits<GraphType>::vertices_size_type VertexIndex;
typedef boost::graph_traits<GraphType>::edges_size_type EdgeIndex;


//edges, weights, source_id, target_id, num_faces
typedef std::tuple<std::vector<std::pair<int, int>>, std::vector<double>, int, int, int> Graph;

// (minX, minY), (maxX, maxY), stepX, stepY, nr_cells
typedef std::tuple<Point, Point, double, double, int> grid;

#endif //PCH_HPP
