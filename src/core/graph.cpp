//
// Created by samuel-berge on 8/11/25.
//
#include "core/graph.hpp"

namespace GRAPH {


    double calculate_length(const Arrangement::Halfedge_const_handle& he) {
        auto p1 = he->source()->point();
        auto p2 = he->target()->point();

        auto dx = p2.x() - p1.x();
        auto dy = p2.y() - p1.y();

        auto squared_length = dx*dx + dy*dy;
        return std::sqrt(CGAL::to_double(squared_length));
    }

    Graph build_graph(const Arrangement& arr, double alpha, Logger &logger) {
        std::vector<std::pair<int, int>> edges;
        std::vector<double> weights;
        int num_faces = arr.number_of_faces()-1; // Exclude the unbounded face
        int source_id = num_faces;
        int target_id = num_faces + 1;

        for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {

            if (fit->is_unbounded() || !fit->has_outer_ccb()) {
                continue; // skip the unbounded face
            }
            double outer_weight = 0.0;
            std::unordered_map<int, double> common_boundaries;

            int face_id = fit->data().id;
            auto curr = fit->outer_ccb();
            do {
                if (!(*curr).twin()->face()->is_unbounded()) {
                    int twin_id= (*curr).twin()->face()->data().id;
                    if (face_id<twin_id) {
                        common_boundaries[twin_id]+=calculate_length(curr);
                    }
                }
                //we don't want orientation of box to influence the solution because it doesn't preserve edge orientation
                //ignore special case with polygon parallel to box because we use offset of box anyway
                else if (!fit->data().belongs_to_poly) {
                    outer_weight = DBL_MAX;
                }
                // in case, face is from poly and lies at outer boundary (very unlikely; has to be parallel to box) and
                // only with outer box with offset 0 possible. then there will be an edge added with w({face, sink})=0
                // and w({face, source})=DBL_MAX so that it is always chosen into the solution
                ++curr;
            } while (curr != fit->outer_ccb());

            for (auto [incident_face, weight] : common_boundaries) {
                edges.emplace_back(face_id, incident_face);
                weights.push_back((1-alpha)*weight);
            }
            if (fit->data().belongs_to_poly) {
                edges.emplace_back(face_id, source_id);
                weights.push_back(DBL_MAX);
                edges.emplace_back(face_id, target_id);
                weights.push_back(0.0);
            }
            else {
                //either outer weight is 0 or infinity, so that outer boundary touching faces are not counted in and
                //the other edges just count in with their area
                edges.emplace_back(face_id, source_id);
                weights.push_back(0);
                edges.emplace_back(face_id, target_id);
                weights.push_back(alpha*(fit->data().area) + (1-alpha)*outer_weight);
            }
        }
        Graph graph (edges, weights, source_id, target_id, num_faces);
        return graph;
    }

}
