#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "pch.hpp"

class algorithm{
public:

  void segments_to_svg(const std::vector<Segment>& segments, const std::string& filename);

  std::pair<std::vector<Segment_w_info>, std::pair<Point, Point>> build_segment_vector(
      const std::vector<SHPLoader::Point2D>& points,
      const std::vector<std::vector<int>>& polygons);

  void add_box(std::vector<Segment_w_info> &segs, std::pair<Point, Point> limits);

  std::vector<Segment> segs_wo_info(const std::vector<Segment_w_info>& segments);

  std::optional<Point> nearest_intersection_in_direction(const Point& origin,
      const Vector& direction,
      const Tree& tree,
      const Segment* self_segment);

  double calculate_length(const Arrangement::Halfedge_const_handle& he);

  std::vector<Segment_w_info> ray_shoot_intersection(const std::vector<Segment_w_info>& segments);



  Arrangement &build_arrangement(std::vector<Segment_w_info> segments);
  
  //return path
  std::tuple<std::vector<std::pair<int, int>>, std::vector<double>, int, int, int> build_graph(const Arrangement& arr, double alpha);

  void AddBidirectionalEdge(GraphType& graph, unsigned int source, unsigned int target, float weight,
                          std::vector<EdgeDescriptor>& reverseEdges, std::vector<float>& capacity);

  const std::vector<Segment> &max_flow(std::tuple<std::vector<std::pair<int, int>>, std::vector<double>, int, int, int> &graph_data, const Arrangement &arr);

  const std::vector<Polygon_2> &output_polygons(const std::vector<int> &groups, VertexDescriptor source, const Arrangement &arr);

  const std::vector<Segment> &output_segs(const std::vector<int> &groups, VertexDescriptor source, const Arrangement &arr);

  void run(std::string input_path);

};

#endif