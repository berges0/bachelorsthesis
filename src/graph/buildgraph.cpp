#include "algorithm.hpp"

std::tuple<std::vector<std::pair<int, int>>, std::vector<double>, int, int, int> algorithm::build_graph(const Arrangement& arr, double alpha){

    std::vector<std::pair<int, int>> edges;
    std::vector<double> weights;
    int nr_edges = 0;
    int nr_vertices = 0;

    // 2. Insert extended segments
    int num_faces = std::distance(arr.faces_begin(), arr.faces_end()) - 1; // Exclude the unbounded face
    int source_id = num_faces;
    int target_id = num_faces + 1; // target is the outer face
    std::cout << "SOURCE , TARGET " << source_id << ", " << target_id << std::endl;
    std::vector<int> vertices(num_faces);
    std::iota(vertices.begin(), vertices.end(), 0);
    Arrangement::Face_const_handle outer_face = arr.unbounded_face();
    int outer_id=outer_face->data().id;
    //std::cout<< "number faces: " << num_faces << std::endl;
    int half_idx=std::ceil(num_faces/2);
    //std::cout<< "half idx: " << half_idx << std::endl;
    int count4=0;
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (fit->is_unbounded()) {
            continue; // skip the unbounded face
        }
        std::vector<std::tuple<int, double>> common_boundaries(0);
        int face_id = fit->data().id;
        if (fit->data().belongs_to_poly) {
            count4++;
        }
        if (fit->has_outer_ccb()) {

            // outer boundary exists
            auto curr = fit->outer_ccb();
            double outer_weight = 0.0;
            do {
                auto twin_face= (*curr).twin()->face();
                if (!twin_face->is_unbounded()) {
                    assert(twin_face->data().id != -1);
                    if (face_id<twin_face->data().id) {
                        bool broke = false;
                        for (int i = 0; i < common_boundaries.size(); ++i) {
                            if (std::get<0>(common_boundaries[i]) == twin_face->data().id) {
                                // storing from perspective of lower id
                                std::get<0>(common_boundaries[i]) = face_id;
                                std::get<1>(common_boundaries[i]) = calculate_length(curr);
                                broke = true;
                                break;
                            }
                        }
                        if(!broke) {
                            common_boundaries.push_back(std::make_tuple(twin_face->data().id, calculate_length(curr)));
                        }
                    }
                }
                else if (!fit->data().belongs_to_poly) {
                    outer_weight += calculate_length(curr);
                }
                ++curr;
            } while (curr != fit->outer_ccb());
            //std::cout << "Face ID: " << face_id << ", the common boundaries is of size " <<  common_boundaries.size()<< std::endl;
            for (auto boundary : common_boundaries) {
                int incident_face_id = std::get<0>(boundary);
                double weight = std::get<1>(boundary);
                    edges.push_back(std::make_pair(face_id, incident_face_id));
                    weights.push_back((1-alpha)*weight);
            }
            if (fit->data().belongs_to_poly) {
                edges.push_back(std::make_pair(face_id, source_id));
                weights.push_back(DBL_MAX);
                edges.push_back(std::make_pair(face_id, target_id));
                weights.push_back(0.0);
            }
            else{
                edges.push_back(std::make_pair(face_id, source_id));
                weights.push_back(0);
                edges.push_back(std::make_pair(face_id, target_id));
                weights.push_back(alpha*(fit->data().area) + (1-alpha)*outer_weight);
            }
        }
    }
    for (int j=0; j< edges.size();j++){
        assert(edges.size()==weights.size());
        //std::cout << "Edge from " << edges[j].first <<" to "<< edges[j].second << " with weight " << weights[j]<<std::endl;

    }
    //std::cout << "Outer face ID: " << outer_id << std::endl;
    return std::make_tuple(edges, weights, source_id, target_id, num_faces);
}



