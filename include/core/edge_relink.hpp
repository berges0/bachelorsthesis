//
// Created by samuel-berge on 8/30/25.
//

#ifndef EDGE_RELINK_HPP
#define EDGE_RELINK_HPP
#include "utils/pch.hpp"


namespace EDGE_RELINK {

    void relink_edges(std::vector<std::vector<Segment_w_info>> &segments);

    std::vector<Point> order_endpoints_along_main_dir(const std::vector<Segment>& segs);

}
#endif //EDGE_RELINK_HPP
