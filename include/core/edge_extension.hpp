//
// Created by samuel-berge on 8/10/25.
//

#ifndef EDGE_EXTENSION_HPP
#define EDGE_EXTENSION_HPP
#include <vector>

#include "utils/pch.hpp"

namespace EDGE_EXTENSION {

    std::vector<Segment> filter_segments(const std::vector<Segment_w_info>& segments);

void add_outer_box(std::vector<Segment_w_info>& segments, double offset);

namespace STANDARD {
const std::vector<Segment_w_info> extension(std::vector<Segment_w_info>& segments);
}

namespace LIMITED{
const std::vector<Segment_w_info> extension(std::vector<Segment_w_info>& segments, double threshold, int th_variant);

const std::vector<Segment_w_info> post_process(const std::vector<Segment_w_info>& segments, const std::vector<bool>& to_prune);
}

    std::optional<Point> first_intersection(const Point& origin,
    const Vector& direction,
    const Tree& tree,
    const Segment* self_segment);

    double get_distance(Point p1, Point p2);

}


#endif //EDGE_EXTENSION_HPP
