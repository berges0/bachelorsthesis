//
// Created by samuel-berge on 8/10/25.
//

#ifndef ARRANGEMENT_HPP
#define ARRANGEMENT_HPP

#include "utils/pch.hpp"

namespace ARRANGEMENT {

    Arrangement build_arrangement(const std::vector<Segment_w_info>& segments);

    void add_edge_data(Arrangement &arr);

    void add_face_data(Arrangement &arr);

    Arrangement build_arrangement_relinked(const std::vector<Segment_w_info>& segments, RTree &rtree, const std::vector<PWH> &polygonswh);

    void add_face_data_relinked(Arrangement &arr, RTree &rtree, const std::vector<PWH> &poly);

    bool test_in_poly(const Point &qp, const RTree &rtree, const std::vector<PWH> &polys);

}

#endif //ARRANGEMENT_HPP
