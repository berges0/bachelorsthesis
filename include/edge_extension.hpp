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
using RTree = bgi::rtree<Segment, bgi::quadratic<16>>;

class EdgeExtension {
public:
    EdgeExtension()=default;
    void buildRTree(const std::vector<SHPLoader::Point2D>& points,
                        const std::vector<std::vector<int>>& polygons);
    std::vector<Segment> queryNearest(double x, double y, int k) const;

private:
    RTree rtree;
};

#endif //EDGE_EXTENSION_HPP
