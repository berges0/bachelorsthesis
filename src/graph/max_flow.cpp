#include "algorithm.hpp"

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

//const std::vector<Polygon_2>
const std::vector<int> &algorithm::max_flow(std::tuple<std::vector<std::pair<int, int>>, std::vector<double>, int, int, int> &graph_data, const Arrangement &arr) {
    GraphType graph;
    int sourceId = std::get<2>(graph_data);
    int sinkId = std::get<3>(graph_data);
      unsigned int numberOfVertices = std::get<4>(graph_data)+2; // +2 for source and sink vertices

      static std::vector<int> groups(numberOfVertices);

      std::vector<EdgeDescriptor> reverseEdges;

      std::vector<float> capacity;

       for (size_t i = 0; i < std::get<0>(graph_data).size(); ++i) {
           unsigned int u_idx = static_cast<unsigned int>(std::get<0>(graph_data)[i].first);
           unsigned int v_idx = static_cast<unsigned int>(std::get<0>(graph_data)[i].second);
           float cap = static_cast<float>(std::get<1>(graph_data)[i]);
           AddBidirectionalEdge(graph, u_idx, v_idx, cap, reverseEdges, capacity);
           AddBidirectionalEdge(graph, v_idx, u_idx, cap, reverseEdges, capacity);
       }



  std::vector<float> residual_capacity(num_edges(graph), 0);

  VertexDescriptor sourceVertex = vertex(sourceId,graph);
  VertexDescriptor sinkVertex = vertex(sinkId,graph);


  std::cout << "Number of vertices " << num_vertices(graph) << std::endl;

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

    int count = 0;
    std::cout << "Source group label " << groups[sourceVertex] << std::endl;
    std::cout << "Sink group label " << groups[sinkVertex] << std::endl;
  for(size_t index=0; index < numberOfVertices - 2; ++index)
  {
    if(groups[index] == groups[sourceVertex])
    {
        count++;
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
    std::cout << "There are "<< count << " faces in the solution" << std::endl;

    std::cout << "Source group " << groups[sourceVertex] << std::endl;
    std::cout << "Sink group " << groups[sinkVertex] << std::endl;
    int counter=0;
    int count3=0;
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (fit->is_unbounded()) continue; // skip unbounded faces
        if (fit->data().belongs_to_poly) {
            count3++;
            std::cout << " IS IN SOLUTION:  "<< (groups[fit->data().id] == groups[sourceVertex]);
            if (!(groups[fit->data().id] == groups[sourceVertex])){
                counter++;
            }
        }
    }
    std::cout << "There are " << counter << " faces out of "<< count3<<" not in the solution although they should be" << std::endl;
    return groups;
}

const std::vector<Polygon_2> &algorithm::output_polygons(const std::vector<int> &groups, const Arrangement &arr) {
    static std::vector<Polygon_2> polygons;
    int sinkId = groups[groups.size() - 1];
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (fit->is_unbounded() || groups[fit->data().id]==0) continue; // skip unbounded faces
        Polygon_2 poly;
        auto circ = fit->outer_ccb();
        auto curr = circ;
        do {
            const auto& source = curr->source()->point();
            const auto& target = curr->target()->point();
            poly.push_back(source);
            ++curr;
        } while (curr != circ);
        if (!poly.is_simple()) {
            std::cerr << "Polygon is not simple!" << std::endl;
            continue;
        }
        polygons.push_back(poly);
    }
    return polygons;
}

const std::vector<Polygon_2> &algorithm::combined_output_polygons(const std::vector<int> &groups, Arrangement &arr) {
    static std::vector<Polygon_2> polygons;
    for (Arrangement::Face_iterator fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (fit->is_unbounded() || !(fit->data().yet_unvisited) || groups[fit->data().id]==0) continue;
            Polygon_2 poly;
            DFS(fit,fit->outer_ccb(), groups,poly);
            polygons.push_back(poly);
    }
    return polygons;
}

void algorithm::DFS(Arrangement::Face_iterator &fit, Arrangement::Halfedge_handle shared_edge ,const std::vector<int> &groups, Polygon_2 &polygon) {
    if (groups[fit->data().id]==0||fit->is_unbounded()||!(fit->data().yet_unvisited)) return;
    fit->data().yet_unvisited=false;
    Arrangement::Ccb_halfedge_circulator circ = fit->outer_ccb();
    Arrangement::Ccb_halfedge_circulator old_start = circ;
    do {
        if (&(*circ) == &(*shared_edge)) {
            // Found the exact starting edge
            break;
        }
        ++circ;
    } while (circ != old_start);
    Arrangement::Ccb_halfedge_circulator new_start = circ;
    do {
        Arrangement::Ccb_halfedge_circulator  twin = circ->twin();
        if (groups[twin->face()->data().id]==0 || twin->face()->is_unbounded()) {
            polygon.push_back(circ->curve().source());
            polygon.push_back(circ->curve().target());
        }
        else if (twin->face()->data().yet_unvisited) {
            Arrangement::Face_iterator fit1 = twin->face();
            DFS(fit1, twin, groups, polygon);
        }
        circ ++;
    }while (new_start != circ);
}
/*
Arrangement::Ccb_halfedge_circulator circ = b->outer_ccb();
Arrangement::Ccb_halfedge_circulator start = circ;
do {
    if (&(*circ) == &(*e->twin())) {
        // Found the exact starting edge
        break;
    }
    ++circ;
} while (circ != start);

// Now `circ` starts at e->twin()
Arrangement::Ccb_halfedge_circulator begin = circ;
do {
    // Traverse b's outer CCB starting from e->twin()
    ++circ;
} while (circ != begin);

*/
const std::vector<Segment> &algorithm::output_segs(const std::vector<int> &groups, VertexDescriptor source, const Arrangement &arr) {
    static std::vector<Segment> segs;
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (fit->is_unbounded()) continue; // skip unbounded faces

        if (groups[fit->data().id]==groups[source]) {
            auto circ = fit->outer_ccb();
            auto curr = circ;
            do {
                const auto& source = curr->source()->point();
                const auto& target = curr->target()->point();
                segs.emplace_back(Segment(source, target));
                ++curr;
            } while (curr != circ);
        }
    }
    return segs;
}
