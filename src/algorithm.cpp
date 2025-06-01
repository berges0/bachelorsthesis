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

    segments_to_svg(segments, "test.svg");

    std::vector<Segment> extended_segments;

    for (const Segment& seg : segments) {
        Line line(seg);
        Point p1 = seg.source();
        Point p2 = seg.target();

        Vector forward = p2 - p1;
        Vector backward = p1 - p2;

        auto hit_forward = nearest_intersection_in_direction(p1, forward, tree, seg);
        auto hit_backward = nearest_intersection_in_direction(p2, backward, tree, seg);

        if (hit_forward) {
            std::cout << "entered" << std::endl;
            extended_segments.emplace_back(p1, *hit_forward);
        }
        if (hit_backward) {
            std::cout << "entered" << std::endl;
            extended_segments.emplace_back(p2, *hit_backward);
        }
    }
    std::cout << tree.size() << std::endl;
    std::cout << extended_segments.size() << std::endl;

    return extended_segments;
}




std::optional<Point> algorithm::nearest_intersection_in_direction(
    const Point& origin,
    const Vector& direction,
    const Tree& tree,
    const Segment& self_segment)
{

    Ray ray(origin, origin + direction);

    auto result = tree.first_intersection(ray);
    if (result) {
        std::cout << "Intersection found\n";
        // unpack intersection as before, return point
    } else {
        std::cout << "No intersection found from origin: "
                  << CGAL::to_double(origin.x()) << ", "
                  << CGAL::to_double(origin.y()) << "\n";
    }
    if (result) {
        if (const Point* ipoint = std::get_if<Point>(&(result->first))) {
            if (*ipoint != self_segment.source() && *ipoint != self_segment.target()) {
                return *ipoint;
            }
        }
    }
    return std::nullopt;
}

// creates a svg, that represents all segments and centers them
void algorithm::segments_to_svg(const std::vector<Segment>& segments, const std::string& filename) {

    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();

    for (const auto& s : segments) {
        min_x = std::min({min_x, CGAL::to_double(s.source().x()), CGAL::to_double(s.target().x())});
        max_x = std::max({max_x, CGAL::to_double(s.source().x()), CGAL::to_double(s.target().x())});
        min_y = std::min({min_y, CGAL::to_double(s.source().y()), CGAL::to_double(s.target().y())});
        max_y = std::max({max_y, CGAL::to_double(s.source().y()), CGAL::to_double(s.target().y())});
    }

    double padding = 10.0; // Adjust based on your coordinate range
    min_x -= padding;
    max_x += padding;
    min_y -= padding;
    max_y += padding;

    double width = max_x - min_x;
    double height = max_y - min_y;

    std::ofstream svg_file(filename);
    svg_file << std::fixed << std::setprecision(10); // HIGH precision output

    svg_file << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
             << "viewBox=\"" << min_x << " " << min_y << " " << width << " " << height << "\" "
             << "width=\"" << 1000 << "\" height=\"" << 1000 * height / width << "\" "
             << "preserveAspectRatio=\"xMidYMid meet\" "
             << "fill=\"none\" stroke=\"black\" stroke-width=\"1\">\n";


    for (const auto& s : segments) {
        svg_file << "<line x1=\"" << s.source().x() << "\" y1=\"" << (max_y -(s.source().y()- min_y))
                 << "\" x2=\"" << s.target().x() << "\" y2=\"" << (max_y - (s.target().y() - min_y))
                 << "\" />\n";
    }

    svg_file << "</svg>\n";

    svg_file.close();
}
