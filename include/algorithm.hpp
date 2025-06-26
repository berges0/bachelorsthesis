#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "pch.hpp"



struct SkipSegment {
  const Segment* to_skip;

  SkipSegment(const Segment* s) : to_skip(s) {}

  bool operator()(const Primitive::Id& id) const {
    const Segment& s = *id;

    // Skip exact same segment (by pointer)
    if (&s == to_skip) return true;

    // Skip segments that share an endpoint with the segment to skip
    return shares_endpoint(s, *to_skip);
  }

  bool shares_endpoint(const Segment& a, const Segment& b) const {
    return a.source() == b.source() || a.source() == b.target() ||
           a.target() == b.source() || a.target() == b.target();
  }
};

struct Segment_w_info{
  Segment seg;
  bool from_poly;
  int poly_id;
  Segment_w_info(const Segment& s, bool fp, int id) : seg(s), from_poly(fp), poly_id(id) {}
};

class algorithm{
public:
  std::vector<Segment_w_info> build_segment_vector(
      const std::vector<SHPLoader::Point2D>& points,
      const std::vector<std::vector<int>>& polygons);

  std::vector<Segment> segs_wo_info(const std::vector<Segment_w_info>& segments);

  std::optional<Point> nearest_intersection_in_direction(const Point& origin,
      const Vector& direction,
      const Tree& tree,
      const Segment* self_segment);

  std::vector<Segment_w_info> ray_shoot_intersection(const std::vector<Segment_w_info>& segments);


  void segments_to_svg(const std::vector<Segment>& segments, const std::string& filename);

  Arrangement &build_arrangement(std::vector<Segment_w_info> segments);

  void run(std::string input_path);

};
#endif