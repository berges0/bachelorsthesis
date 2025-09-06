//
// Created by samuel-berge on 8/30/25.
//

#ifndef SUBSTITUTE_EDGES_HPP
#define SUBSTITUTE_EDGES_HPP
#include "utils/pch.hpp"


namespace SUBSTITUTE_EDGES {

    void relink_edges(std::vector<std::vector<Segment_w_info>> &segments);

    std::vector<Point> order_endpoints_along_main_dir(const std::vector<Segment>& segs);

    RTree build_r_tree(const std::vector<PWH>& polygons);

    void connect_outer_points(std::vector<std::vector<Segment_w_info>> &segments);
}
#endif //SUBSTITUTE_EDGES_HPP
