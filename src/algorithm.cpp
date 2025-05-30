#include "algorithm.hpp"

const std::vector<Segment> &algorithm::build_segment_vector(const std::vector<SHPLoader::Point2D>& points,
                    const std::vector<std::vector<int>>& polygons) {
    std::vector <Segment> segments;
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

std::optional<Point> algorithm::first_intersection(const Ray& ray, const Tree& tree, const Segment& self) {
    auto result = tree.first_intersection(ray);

    if (!result) return std::nullopt;
    const auto& intersection_variant = result->first; // boost::variant<Point, Segment>

    if (const Point* p = boost::get<Point>(&intersection_variant)) {
        return *p;
    } else if (const Segment* s = boost::get<Segment>(&intersection_variant)) {
        Point src = ray.source();
        double d1 = CGAL::squared_distance(src, s->source());
        double d2 = CGAL::squared_distance(src, s->target());
        return d1 < d2 ? s->source() : s->target();
    }
    return std::nullopt;
}

std::vector<Segment> algorithm::extend_segments(const std::vector<Segment>& segments) {
    Tree tree(segments.begin(), segments.end());
    tree.accelerate_distance_queries();

    std::vector<Segment> extended;

    for (const Segment& seg : segments) {
        Point p1 = seg.source();
        Point p2 = seg.target();

        Vector dir = p2 - p1;
        Vector rev = -dir;

        Ray forward(p2, dir);
        Ray backward(p1, rev);

        auto f_hit = first_intersection(forward, tree, seg);
        auto b_hit = first_intersection(backward, tree, seg);

        Point new_start = b_hit.value_or(p1);
        Point new_end   = f_hit.value_or(p2);

        extended.emplace_back(new_start, new_end);
    }

    return extended;
}

