#include "algorithm.hpp"

std::tuple<std::vector<std::pair<int, int>>, std::vector<double>, int, int, int> algorithm::build_graph(const Arrangement& arr, double alpha){

    std::vector<std::pair<int, int>> edges;
    std::vector<double> weights;   
    int nr_edges = 0;
    int nr_vertices = 0;

    // 2. Insert extended segments
    int num_faces = std::distance(arr.faces_begin(), arr.faces_end());
    int source_id = num_faces + 1;
    int target_id = num_faces + 2; // target is the outer face
    std::cout << "SOURCE , TARGET " << source_id << ", " << target_id << std::endl;
    std::vector<int> vertices(num_faces);
    std::iota(vertices.begin(), vertices.end(), 0);
    Arrangement::Face_const_handle outer_face = arr.unbounded_face();
    int outer_id=outer_face->data().id;
    //std::cout<< "number faces: " << num_faces << std::endl;
    int half_idx=std::ceil(num_faces/2);
    //std::cout<< "half idx: " << half_idx << std::endl;

    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if(fit->data().id == outer_id) {
            continue; // skip the outer face
        }
        std::vector<std::tuple<int, double>> common_boundaries(0);
        int face_id = fit->data().id;
        if (fit->has_outer_ccb()) {
            // outer boundary exists
            auto curr = fit->outer_ccb();
            double outer_weight = 0.0;
            do {
                int incident_face_id = (*curr).twin()->face()->data().id ;
                if (incident_face_id != outer_id) {
                    if (face_id<=half_idx){
                        bool broke = false;
                        for (int i = 0; i < common_boundaries.size(); ++i) {
                            if (std::get<0>(common_boundaries[i]) == incident_face_id) {
                                // storing from perspective of lower id
                                std::get<0>(common_boundaries[i]) = face_id;
                                std::get<1>(common_boundaries[i]) = calculate_length(curr);
                                broke = true;
                                break;
                            }
                        }
                        if(!broke) {
                            common_boundaries.push_back(std::make_tuple(incident_face_id, calculate_length(curr)));
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
    return std::make_tuple(edges, weights, source_id, target_id, num_faces+2);
}



void algorithm::AddBidirectionalEdge(GraphType& graph, unsigned int source, unsigned int target, float weight,
                          std::vector<EdgeDescriptor>& reverseEdges, std::vector<float>& capacity)
{
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

void algorithm::max_flow(std::tuple<std::vector<std::pair<int, int>>, std::vector<double>, int, int, int> &graph_data) {
    GraphType graph;

      unsigned int numberOfVertices = std::get<4>(graph_data);
      std::vector<int> groups(numberOfVertices);

      std::vector<EdgeDescriptor> reverseEdges;

      std::vector<float> capacity;

       for (size_t i = 0; i < std::get<0>(graph_data).size(); ++i) {
           unsigned int u_idx = static_cast<unsigned int>(std::get<0>(graph_data)[i].first);
           unsigned int v_idx = static_cast<unsigned int>(std::get<0>(graph_data)[i].second);
           float cap = static_cast<float>(std::get<1>(graph_data)[i]);
           AddBidirectionalEdge(graph, u_idx, v_idx, cap, reverseEdges, capacity);
           AddBidirectionalEdge(graph, v_idx, u_idx, cap, reverseEdges, capacity);
       }


  int sourceId = 9;
  int sinkId = 10;

  std::vector<float> residual_capacity(num_edges(graph), 0);

  VertexDescriptor sourceVertex = vertex(std::get<4>(graph_data)-2,graph);
  VertexDescriptor sinkVertex = vertex(std::get<4>(graph_data)-1,graph);

  // There should be 3*3 + 2 = 11 nodes
  std::cout << "Number of vertices " << num_vertices(graph) << std::endl;

  // There should be 2*(12) + 2*2*9 = 60 edges
  std::cout << "Number of edges " << num_edges(graph) << std::endl;

  boost::boykov_kolmogorov_max_flow(graph,
      boost::make_iterator_property_map(&capacity[0], get(boost::edge_index, graph)),
      boost::make_iterator_property_map(&residual_capacity[0], get(boost::edge_index, graph)),
      boost::make_iterator_property_map(&reverseEdges[0], get(boost::edge_index, graph)),
      boost::make_iterator_property_map(&groups[0], get(boost::vertex_index, graph)),
      get(boost::vertex_index, graph),
      sourceVertex,
      sinkVertex);

  // Display the segmentation
//  for(size_t index=0; index < groups.size(); ++index)
//  {
//       std::cout << "Vertex " << index << " is in group " << groups[index] << std::endl;
//  }

  std::cout << "Source group " << groups[sourceVertex] << std::endl;
  std::cout << "Sink group " << groups[sinkVertex] << std::endl;

  for(size_t index=0; index < numberOfVertices - 2; ++index)
  {
    if(groups[index] == groups[sourceVertex])
    {
        std::cout << "Vertex " << index << " is attached to the source. (group " << groups[index] << ")" << std::endl;
    }
    else if(groups[index] == groups[sinkVertex])
    {
        std::cout << "Vertex " << index << " is attached to the sink. (group " << groups[index] << ")" << std::endl;
    }
    else
    {
        std::cerr << "Vertex " << index << " is not attached to the source or sink! (group " <<  groups[index] << ")" << std::endl;
    }
  }

}

/*

typedef boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS,
    boost::property<boost::vertex_predecessor_t, int>,                                           // vertex properties
    boost::property<boost::edge_capacity_t, int,                 // edge capacity
        boost::property<boost::edge_residual_capacity_t, int,    // residual capacity
            boost::property<boost::edge_reverse_t,
                           boost::adjacency_list_traits<
                               boost::vecS, boost::vecS, boost::directedS
                           >::edge_descriptor>>>> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;


void algorithm::max_flow(std::tuple<std::vector<std::pair<int, int>>, std::vector<double>, int, int, int> &graph_data) {
    Graph g;
    auto capacity = boost::get(boost::edge_capacity, g);
    auto rev = boost::get(boost::edge_reverse, g);
    Vertex src_dscrpt, tar_dscrpt;
    std::vector<Vertex> vertices(std::get<4>(graph_data));
    for (int i = 0; i < std::get<4>(graph_data); ++i) {
        vertices[i] = boost::add_vertex(g);
        if (i == std::get<4>(graph_data)-2) {
            src_dscrpt = vertices.back();
        }
        if (i == std::get<4>(graph_data)-1) {
            tar_dscrpt = vertices.back();
        }
    }
    for (size_t i = 0; i < std::get<0>(graph_data).size(); ++i) {
        int u_idx = std::get<0>(graph_data)[i].first;
        int v_idx = std::get<0>(graph_data)[i].second;
        double cap = std::get<1>(graph_data)[i];

        bool success;

        // Edge for u→v and its reverse
        Edge e_uv, e_uv_rev;
        boost::tie(e_uv, success) = boost::add_edge(vertices[u_idx], vertices[v_idx], g);
        boost::tie(e_uv_rev, success) = boost::add_edge(vertices[v_idx], vertices[u_idx], g);

        capacity[e_uv] = cap;
        capacity[e_uv_rev] = 0;
        rev[e_uv] = e_uv_rev;
        rev[e_uv_rev] = e_uv;

        // Edge for v→u and its reverse (if you want undirected behavior with capacity both ways)
        Edge e_vu, e_vu_rev;
        boost::tie(e_vu, success) = boost::add_edge(vertices[v_idx], vertices[u_idx], g);
        boost::tie(e_vu_rev, success) = boost::add_edge(vertices[u_idx], vertices[v_idx], g);

        capacity[e_vu] = cap;
        capacity[e_vu_rev] = 0;
        rev[e_vu] = e_vu_rev;
        rev[e_vu_rev] = e_vu;

    }
    std::cout << "Graph has " << boost::num_vertices(g) << " vertices and "
              << boost::num_edges(g) << " edges." << std::endl;

    auto maxflow = boost::boykov_kolmogorov_max_flow(g, src_dscrpt, tar_dscrpt);
}
*/