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

const std::vector<Polygon_with_holes_2> algorithm::combined_output_polygons(const std::vector<int> &groups,
                                                                                     Arrangement &arr) {
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (fit->is_unbounded()) {
            continue; // skip unbounded faces
        }
        Polygon_2 test_poly;
        auto current = fit->outer_ccb();
        do {
            test_poly.push_back(current->source()->point());
            current++;
        }while(current != fit->outer_ccb());
        std::cout << "Face is counterclockwise " << test_poly.is_counterclockwise_oriented() << std::endl;
    }



    std::vector<Polygon_with_holes_2> polygons;
    std::vector<Polygon_2> outer_boundaries;
    std::vector<Polygon_2> holes;
    for (Arrangement::Edge_iterator eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
        if (!eit->data().visited && groups[eit->face()->data().id] != 0 && groups[eit->twin()->face()->data().id] == 0) {
            Arrangement::Halfedge_handle eit_handle = eit;
            Polygon_2 poly = merge_contiguous_faces(eit_handle, groups, arr);
            assert(poly.is_simple());
            if (poly.is_counterclockwise_oriented()) {
                outer_boundaries.emplace_back(poly);
            } else {
                holes.emplace_back(poly);
            }
        }
    }
    std::vector<CGAL::Bbox_2> outer_bboxes;
    for (const auto& poly : outer_boundaries) {
        outer_bboxes.emplace_back(poly.bbox());
    }

    std::map<int, std::vector<Polygon_2>> outer_to_holes;
    for (int i = 0; i < holes.size(); ++i) {
        const Polygon_2& hole = holes[i];
        Point test_point = *hole.vertices_begin();
        CGAL::Bbox_2 point_bbox = test_point.bbox();

        for (int j = 0; j < outer_boundaries.size(); ++j) {
            if (!CGAL::do_overlap(point_bbox, outer_bboxes[j]))
                continue;

            CGAL::Bounded_side side = CGAL::bounded_side_2(
                outer_boundaries[j].vertices_begin(),
                outer_boundaries[j].vertices_end(),
                test_point,
                Kernel());

            if (side == CGAL::ON_BOUNDED_SIDE) {
                outer_to_holes[j].emplace_back(holes[i]);  // hole i is inside outer j
                break;
            }
        }
    }

    std::vector<bool> already_added (outer_boundaries.size(), false);
    for (const auto &outer_pair : outer_to_holes) {
        int compare = polygons.size();
        Polygon_with_holes_2 pwh(outer_boundaries[outer_pair.first], outer_pair.second.begin(), outer_pair.second.end());
        polygons.emplace_back(pwh);
        already_added[outer_pair.first] = true;
    }
    for (int outer_boundary = 0; outer_boundary < outer_boundaries.size(); ++outer_boundary) {
        if (!already_added[outer_boundary]) {
            int compare = polygons.size();

            Polygon_with_holes_2 pwh(outer_boundaries[outer_boundary]);
            polygons.emplace_back(pwh);

        }
    }
    return polygons;
}

int glob_count = 0;

const Polygon_2 algorithm::merge_contiguous_faces(Arrangement::Halfedge_handle &edge, const std::vector<int> &groups, Arrangement &arr ) {
    int count = 0;
    assert(groups[edge->face()->data().id] != 0 && groups[edge->twin()->face()->data().id] == 0);
    Polygon_2 polygon;
    do {

        if (groups[edge->twin()->face()->data().id] == 0 || edge->twin()->face()->is_unbounded()) {
            //assert(edge->data().visited==false);
            polygon.push_back(edge->source()->point());
            count++;
            if (edge->target()->point() == *(polygon.vertices_begin())) {
                polygon.push_back(edge->target()->point());
                return polygon; // finished

            }
            edge = edge->next();
        }
        else if (groups[edge->twin()->face()->data().id] != 0) {
            edge = edge->twin()->next();
        }
    }while (true);

}


/*
const std::vector<Polygon_2> algorithm::get_holes(Arrangement::Halfedge_handle &edge, const std::vector<int> &groups, Arrangement &arr ) {
    Arrangement::Halfedge_handle current = edge;
    Arrangement::Halfedge_handle start = current;

    do {

    }while

    do {
        while (current->data().outer == true && current->data().visited == true) {
            current = current->next();
        }
        auto start = current->prev();
        do{
            if (groups[current->face()->data().id]==0) {
                
            }
            current = current->next();
        }while (current->prev() != start);


    }
}*/


const std::vector<Polygon_2> &algorithm::output_polygons(const std::vector<int> &groups, const Arrangement &arr) {
    static std::vector<Polygon_2> polygons;
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
