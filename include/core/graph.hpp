//
// Created by samuel-berge on 8/10/25.
//

#ifndef GRAPH_HPP
#define GRAPH_HPP
#include <vector>
#include <bits/stl_pair.h>

#include "utils/pch.hpp"

namespace GRAPH {

    Graph &build_graph(const Arrangement& arr, double alpha);

    void update_common_boundaries(int face_1, int face_2, double weight, std::vector<std::pair<int, double>> &common_boundaries);

    void add_polyface(std::vector<std::pair<int, int>> &edges, std::vector<double> &weights, int s, int t, int face_id);

    void add_non_polyface(std::vector<std::pair<int, int>> &edges, std::vector<double> &weights, int s, int t, int face_id,
    double alpha, double outer_weight, double area);


}
#endif //GRAPH_HPP
