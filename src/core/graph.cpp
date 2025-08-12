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

    Graph &build_graph(const Arrangement& arr, double alpha) {
        std::vector<std::pair<int, int>> edges;
        std::vector<double> weights;

        int num_faces = std::distance(arr.faces_begin(), arr.faces_end()) - 1; // Exclude the unbounded face
        int source_id = num_faces;
        int target_id = num_faces + 1;

        for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
            if (fit->is_unbounded()) {
                continue; // skip the unbounded face
            }
            std::vector<std::pair<int, double>> common_boundaries(0);
            int face_id = fit->data().id;

            if (fit->has_outer_ccb()) {
                // outer boundary exists
                auto curr = fit->outer_ccb();
                double outer_weight = 0.0;
                do {
                    auto twin_face= (*curr).twin()->face();
                    if (!twin_face->is_unbounded()) {
                        int twin_id= twin_face->data().id;
                        update_common_boundaries(face_id,twin_id, calculate_length(curr), common_boundaries);
                    }
                    else if (!fit->data().belongs_to_poly) {
                        outer_weight += calculate_length(curr);
                    }
                    ++curr;
                } while (curr != fit->outer_ccb());

                for (auto boundary : common_boundaries) {
                    int incident_face_id = std::get<0>(boundary);
                    double weight = std::get<1>(boundary);
                    edges.push_back(std::make_pair(face_id, incident_face_id));
                    weights.push_back((1-alpha)*weight);
                }

                if (fit->data().belongs_to_poly) {
                    add_polyface(edges, weights, source_id, target_id, face_id);
                }
                else{
                    add_non_polyface(edges, weights, source_id, target_id, face_id,
                        alpha, outer_weight, (fit->data().area));
                }
            }
        }

        static Graph graph(edges, weights, source_id, target_id, num_faces);
        return graph;
    }

    void update_common_boundaries(int face_1, int face_2, double weight, std::vector<std::pair<int, double>> &common_boundaries) {
        assert(face_1 != -1 && face_2 != -1);
        if (face_1<face_2) {
            bool broke = false;
            for (int i = 0; i < common_boundaries.size(); ++i) {
                if (std::get<0>(common_boundaries[i]) == face_2) {
                    // storing from perspective of lower id
                    std::get<0>(common_boundaries[i]) = face_1;
                    std::get<1>(common_boundaries[i]) += weight;
                    broke = true;
                    break;
                }
            }
            if(!broke) {
                common_boundaries.push_back(std::pair(face_2, weight));
            }
        }
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
