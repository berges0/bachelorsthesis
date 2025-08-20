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

    Graph build_graph(const Arrangement& arr, double alpha) {
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
                else if (!fit->data().belongs_to_poly) {
                    outer_weight += calculate_length(curr);
                }
                ++curr;
            } while (curr != fit->outer_ccb());

            for (auto [incident_face, weight] : common_boundaries) {
                edges.push_back(std::make_pair(face_id, incident_face));
                weights.push_back((1-alpha)*weight);
            }
            if (fit->data().belongs_to_poly) {
                add_polyface(edges, weights, source_id, target_id, face_id);
            }
            else {
                add_non_polyface(edges, weights, source_id, target_id, face_id,
                    alpha, outer_weight, (fit->data().area));
            }
        }
        Graph graph (edges, weights, source_id, target_id, num_faces);
        return graph;
    }


    void add_polyface(std::vector<std::pair<int, int>> &edges, std::vector<double> &weights, int s, int t, int face_id) {
        edges.push_back(std::make_pair(face_id, s));
        weights.push_back(DBL_MAX);
        edges.push_back(std::make_pair(face_id, t));
        weights.push_back(0.0);
    }

    void add_non_polyface(std::vector<std::pair<int, int>> &edges, std::vector<double> &weights, int s, int t, int face_id,
        double alpha, double outer_weight, double area) {
        edges.push_back(std::make_pair(face_id, s));
        weights.push_back(0);
        edges.push_back(std::make_pair(face_id, t));
        weights.push_back(alpha*area + (1-alpha)*outer_weight);
    }
}
