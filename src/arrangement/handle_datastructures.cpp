#include "algorithm.hpp"

std::pair<std::vector<Segment_w_info>, std::pair<Point, Point>> algorithm::build_segment_vector(const std::vector<SHPLoader::Point2D>& points,
                    const std::vector<std::vector<int>>& polygons) {
    std::vector<Segment_w_info> segments;
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    std::cout <<"There are "<<polygons.size() <<"polygons..." << std::endl;

    int id=0;
    for (const auto &poly:polygons) {
        int n = poly.size();
        if (n < 2) continue;

        for (int i = 0; i < n; ++i) {
            const auto& p1 = points[poly[i]];
            const auto& p2 = points[poly[(i + 1) % n]]; // Polygon schließen
            if (p1.x < minX) minX = p1.x;
            if (p1.x > maxX) maxX = p1.x;
            if (p1.y < minY) minY = p1.y;
            if (p1.y > maxY) maxY = p1.y;
            Point pt1(p1.x, p1.y);
            Point pt2(p2.x, p2.y);
            Segment_w_info edge(Segment(pt1, pt2), true, id);
            segments.push_back(edge);
        }
        ++id;
    }
    std::pair<Point, Point> limits = {Point(minX, minY), Point(maxX, maxY)};
    std::vector<Curve> curves;
    Arrangement arr;
    for (const auto& seg : segments) {
        Curve c(seg.seg, seg.from_poly);
        //std::cout << "Added segment from" << seg.from_poly << " polygon " << std::endl;
        curves.push_back(c);
    }

    CGAL::insert(arr, curves.begin(), curves.end());
    int num_faces = std::distance(arr.faces_begin(), arr.faces_end()) - 1; // Exclude the unbounded face
    int should_be = polygons.size();

    return {segments, limits};
}

std::vector<Segment> algorithm::segs_wo_info(const std::vector<Segment_w_info>& segments) {
    std::vector<Segment> just_segments;
    for (const auto& seg : segments) {
        just_segments.emplace_back(seg.seg);
    }
    return just_segments;
}
