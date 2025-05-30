//
// Created by berges0 on 5/29/25.
//

#ifndef EDGE_EXTENSION_HPP
#define EDGE_EXTENSION_HPP

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <vector>
#include "instance_loader/BasicDataStructures.hpp"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using Point = bg::model::d2::point_xy<double>;
using Segment = bg::model::segment<Point>;


struct MySegment {
  Segment geometry;
  std::vector<std::size_t> polygon_ids; // or just one ID if only one polygon

  MySegment(Point const& p1, Point const& p2, std::vector<std::size_t> ids)
      : geometry(p1, p2), polygon_ids(std::move(ids)) {}
};
namespace boost { namespace geometry { namespace traits {

template<> struct tag<MySegment> { using type = segment_tag; };
template<> struct point_type<MySegment> { using type = Point; };

template<std::size_t Index>
struct indexed_access<MySegment, Index, 0> {
    static double get(MySegment const& s) {
        return boost::geometry::get<Index>(s.geometry.first);
    }
    static void set(MySegment& s, double value) {
        boost::geometry::set<Index>(s.geometry.first, value);
    }
};

template<size_t Index>
struct indexed_access<MySegment, Index, 1> {
    static double get(MySegment const& s) {
        return boost::geometry::get<Index>(s.geometry.second);
    }
    static void set(MySegment& s, double value) {
        boost::geometry::set<Index>(s.geometry.second, value);
    }
};

}}}

struct SegmentTranslator {
    using result_type = boost::geometry::model::segment<Point> const&;

    result_type operator()(MySegment const& s) const {
        return s.geometry;
    }
};

using RTree = bgi::rtree<MySegment, bgi::quadratic<16>, SegmentTranslator>;

class EdgeExtension {
public:
    EdgeExtension()=default;
    void buildRTree(const std::vector<SHPLoader::Point2D>& points,
                        const std::vector<std::vector<int>>& polygons);
    std::vector<MySegment> queryNearest(double x, double y, int k) const;

private:
    RTree rtree;
};

#endif //EDGE_EXTENSION_HPP
