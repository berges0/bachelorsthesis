#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/intersections.h>
#include <vector>
#include <iostream>
#include <optional>
#include "instance_loader/BasicDataStructures.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel; // 2D kernel
using Point        = Kernel::Point_2;
using Segment      = Kernel::Segment_2;
using Ray          = Kernel::Ray_2;
using Vector       = Kernel::Vector_2;

using SegmentIterator = std::vector<Segment>::const_iterator;
using Primitive       = CGAL::AABB_segment_primitive<Kernel, SegmentIterator>;
using Traits          = CGAL::AABB_traits<Kernel, Primitive>;
using Tree            = CGAL::AABB_tree<Traits>;

class algorithm{
    public:
    algorithm() = default;

    const std::vector<Segment> &build_segment_vector(const std::vector<SHPLoader::Point2D>& points,
                        const std::vector<std::vector<int>>& polygons);

    // Computes first intersection of a ray with other segments in the tree
    std::optional<Point> first_intersection(const Ray& ray, const Tree& tree, const Segment& self);

    // Extends each segment in both directions until it hits another segment
    std::vector<Segment> extend_segments(const std::vector<Segment>& segments);
};
#endif //ALGORITHM_HPP
