//
// Created by samuel-berge on 8/10/25.
//

#ifndef MAX_FLOW_HPP
#define MAX_FLOW_HPP

#include "utils/pch.hpp"

namespace MAX_FLOW {

    const std::vector<bool> &max_flow(Graph &graph_data, const Arrangement &arr);

    void add_bidirectional_edge(GraphType& graph, unsigned int source, unsigned int target, float weight,
          std::vector<EdgeDescriptor>& reverseEdges, std::vector<float>& capacity);

    const std::vector<bool> &solution_vector(const std::vector<int> &groups);

}


#endif //MAX_FLOW_HPP
