//
// Created by samuel-berge on 8/10/25.
//

#ifndef ARRANGEMENT_HPP
#define ARRANGEMENT_HPP

#include "utils/pch.hpp"
#include "utils/logger.hpp"

namespace ARRANGEMENT {

    Arrangement build_arrangement(const std::vector<Segment_w_info>& segments, Logger &logger);

    void add_edge_data(Arrangement &arr, Logger &logger);

    void add_face_data(Arrangement &arr, Logger &logger);

    Arrangement build_arrangement_relinked(const std::vector<Segment_w_info>& segments, RTree &rtree, const std::vector<PWH> &polygonswh,
        Logger &logger);

    void add_face_data_relinked(Arrangement &arr, RTree &rtree, const std::vector<PWH> &poly, Logger &logger);

    bool test_in_poly(const Point &qp, const RTree &rtree, const std::vector<PWH> &polys);

}

#endif //ARRANGEMENT_HPP
