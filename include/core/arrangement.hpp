//
// Created by samuel-berge on 8/10/25.
//

#ifndef ARRANGEMENT_HPP
#define ARRANGEMENT_HPP

#include "utils/pch.hpp"

namespace ARRANGEMENT {

    Arrangement build_arrangement(std::vector<Segment_w_info>& segments);

    void add_edge_data(Arrangement &arr);

    void add_face_data(Arrangement &arr);


}

#endif //ARRANGEMENT_HPP
