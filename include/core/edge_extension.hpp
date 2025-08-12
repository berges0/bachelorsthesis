//
// Created by samuel-berge on 8/10/25.
//

#ifndef EDGE_EXTENSION_HPP
#define EDGE_EXTENSION_HPP
#include <vector>

#include "utils/pch.hpp"

namespace EDGE_EXTENSION {

    std::vector<Segment> filter_segments(const std::vector<Segment_w_info>& segments);

    void add_outer_box(std::vector<Segment_w_info>& segments);

    const std::vector<Segment_w_info> &edge_extension(std::vector<Segment_w_info>& segments, std::string version);

    namespace STANDARD {
        const std::vector<Segment_w_info> &extension(std::vector<Segment_w_info>& segments);
    }

    namespace LIMITED{
        const std::vector<Segment_w_info> &extension(std::vector<Segment_w_info>& segments);

        void post_process(std::vector<Segment_w_info>& segments);
    }

    std::optional<Point> first_intersection(const Point& origin,
    const Vector& direction,
    const Tree& tree,
    const Segment* self_segment);
}


#endif //EDGE_EXTENSION_HPP
