//
// Created by samuel-berge on 8/11/25.
//
#include "core/max_flow.hpp"


namespace MAX_FLOW {

    const std::vector<bool> max_flow(Graph &graph_data, const Arrangement &arr) {
        GraphType graph;
        int sourceId = std::get<2>(graph_data);
        int sinkId = std::get<3>(graph_data);
        unsigned int numberOfVertices = std::get<4>(graph_data)+2; // +2 for source and sink vertices

        std::vector<int> groups(numberOfVertices);

        std::vector<EdgeDescriptor> reverseEdges;

        std::vector<float> capacity;

        for (size_t i = 0; i < std::get<0>(graph_data).size(); ++i) {
            unsigned int u_idx = static_cast<unsigned int>(std::get<0>(graph_data)[i].first);
            unsigned int v_idx = static_cast<unsigned int>(std::get<0>(graph_data)[i].second);
            float cap = static_cast<float>(std::get<1>(graph_data)[i]);
            add_bidirectional_edge(graph, u_idx, v_idx, cap, reverseEdges, capacity);
            add_bidirectional_edge(graph, v_idx, u_idx, cap, reverseEdges, capacity);
        }

        std::vector<float> residual_capacity(num_edges(graph), 0);

        VertexDescriptor sourceVertex = vertex(sourceId,graph);
        VertexDescriptor sinkVertex = vertex(sinkId,graph);

        boost::boykov_kolmogorov_max_flow(graph,
            boost::make_iterator_property_map(&capacity[0], get(boost::edge_index, graph)),
            boost::make_iterator_property_map(&residual_capacity[0], get(boost::edge_index, graph)),
            boost::make_iterator_property_map(&reverseEdges[0], get(boost::edge_index, graph)),
            boost::make_iterator_property_map(&groups[0], get(boost::vertex_index, graph)),
            get(boost::vertex_index, graph),
            sourceVertex,
            sinkVertex);

        return solution_vector(groups);
    }

    void add_bidirectional_edge(GraphType& graph, unsigned int source, unsigned int target, float weight,
        std::vector<EdgeDescriptor>& reverseEdges, std::vector<float>& capacity) {
            // Add edges between grid vertices. We have to create the edge and the reverse edge,
            // then add the reverseEdge as the corresponding reverse edge to 'edge', and then add 'edge'
            // as the corresponding reverse edge to 'reverseEdge'
            int nextEdgeId = num_edges(graph);

            EdgeDescriptor edge;
            bool inserted;

            boost::tie(edge,inserted) = add_edge(source, target, nextEdgeId, graph);
            if(!inserted)
            {
                std::cerr << "Not inserted!" << std::endl;
            }
            EdgeDescriptor reverseEdge = add_edge(target, source, nextEdgeId + 1, graph).first;
            reverseEdges.push_back(reverseEdge);
            reverseEdges.push_back(edge);
            capacity.push_back(weight);
            capacity.push_back(weight);
    }

    const std::vector<bool> solution_vector(const std::vector<int> &groups){
        std::vector<bool> solution(groups.size(), false);
        for(size_t index=0; index < groups.size() - 2; ++index) {
            if(groups[index] != 0) {
                solution[index] = true;
            }
        }
        return solution;
    }

}
