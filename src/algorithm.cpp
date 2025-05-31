#include "algorithm.hpp"

std::vector<Segment> algorithm::build_segment_vector(const std::vector<SHPLoader::Point2D>& points,
                    const std::vector<std::vector<int>>& polygons) {
    std::vector<Segment> segments;
    for (const auto &poly:polygons) {
        int n = poly.size();
        if (n < 2) continue;

        for (int i = 0; i < n; ++i) {
            const auto& p1 = points[poly[i]];
            const auto& p2 = points[poly[(i + 1) % n]]; // Polygon schließen

            Point pt1(p1.x, p1.y);
            Point pt2(p2.x, p2.y);
            Segment edge(pt1, pt2);
            segments.push_back(edge);
        }
    }
    return segments;
}

std::vector<Segment> algorithm::ray_shoot_intersection(std::vector<Segment> segments) {
    Tree tree(segments.begin(), segments.end());
    tree.build();
    std::cout << tree.size() << std::endl;
    return segments;
}







