#include "algorithm.hpp"
#include <cmath>

std::vector<Segment_w_info> algorithm::build_segment_vector(const std::vector<SHPLoader::Point2D>& points,
                    const std::vector<std::vector<int>>& polygons) {
    std::vector<Segment_w_info> segmentss;
    int id=0;
    for (const auto &poly:polygons) {
        int n = poly.size();
        if (n < 2) continue;

        for (int i = 0; i < n; ++i) {
            const auto& p1 = points[poly[i]];
            const auto& p2 = points[poly[(i + 1) % n]]; // Polygon schließen

            Point pt1(p1.x, p1.y);
            Point pt2(p2.x, p2.y);
            Segment_w_info edge(Segment(pt1, pt2), true, id);
            segmentss.push_back(edge);
        }
        ++id;
    }
    return segmentss;
}

std::vector<Segment> algorithm::segs_wo_info(const std::vector<Segment_w_info>& segments) {
    std::vector<Segment> just_segments;
    for (const auto& seg : segments) {
        just_segments.emplace_back(seg.seg);
    }
    return just_segments;
}

std::optional<Point> algorithm::nearest_intersection_in_direction(const Point& origin,
    const Vector& direction,
    const Tree& tree,
    const Segment* self_segment){
    Ray ray(origin, origin + direction);
    SkipSegment skip(self_segment);

    auto result = tree.first_intersection(ray,skip);
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
            if (*ipoint != self_segment->source() && *ipoint != self_segment->target()) {
                std::cout<<"and it is not itself but "<<origin<< ", "<< *ipoint <<std::endl;
                return *ipoint;
            }
        }
    }
    return std::nullopt;
}

std::vector<Segment_w_info> algorithm::ray_shoot_intersection(const std::vector<Segment_w_info>& segments) {
    std::vector<Segment> just_segments = segs_wo_info(segments);
    Tree tree(just_segments.begin(), just_segments.end());
    tree.build();

    segments_to_svg(segs_wo_info(segments), "test.svg");

    std::vector<Segment_w_info> all_segments = segments;

    for (const Segment& segment : just_segments) {
        Line line(segment);
        Point p1 = segment.source();
        Point p2 = segment.target();
        std::cout << p1.x() << ", "<< p1.y() << " / " << p2.x() <<", "<< p2.y() << "\n";

        Vector forward = p1 - p2;
        Vector backward = p2 - p1;

        auto hit_forward = nearest_intersection_in_direction(p1, forward, tree, &segment);
        auto hit_backward = nearest_intersection_in_direction(p2, backward, tree, &segment);

        if (hit_forward) {
            all_segments.emplace_back(Segment_w_info(Segment(p1, *hit_forward), false, -1));
            std::cout<< "EMPLACED " << all_segments.back().seg.source().x() << " " << all_segments.back().seg.source().y()<<", "<< all_segments.back().seg.target().x()<<" "<<all_segments.back().seg.target().y()<<std::endl;
        }
        if (hit_backward) {
            all_segments.emplace_back(Segment_w_info(Segment(p2, *hit_backward),false, -1));
            std::cout<< "EMPLACED " << all_segments.back().seg.source().x() << " " << all_segments.back().seg.source().y()<<", "<< all_segments.back().seg.target().x()<<" "<<all_segments.back().seg.target().y()<<std::endl;
        }
    }
    return all_segments;
}


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
    double pixel_width = 1000.0;
    double aspect_ratio = height / width;
    double pixel_height = pixel_width * aspect_ratio;
    double scale = pixel_width / width; // pixel per unit
    double stroke_width_units = 1.0 / scale; // 1 screen pixel = X world units

    svg_file << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
             << "viewBox=\"" << min_x << " " << min_y << " " << width << " " << height << "\" "
             << "width=\"" << pixel_width << "\" height=\"" << pixel_height << "\" "
             << "preserveAspectRatio=\"xMidYMid meet\" "
             << "fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width_units << "\">\n";


    for (const auto& s : segments) {
        svg_file << "<line x1=\"" << s.source().x() << "\" y1=\"" << (max_y -(s.source().y()- min_y))
                 << "\" x2=\"" << s.target().x() << "\" y2=\"" << (max_y - (s.target().y() - min_y))
                 << "\" />\n";
    }

    svg_file << "</svg>\n";

    svg_file.close();
}

Arrangement algorithm::build_arrangement(std::vector<Segment_w_info> segments) {
    Arrangement arr;
    for (const Segment_w_info& seg : segments) {
        CGAL::insert(arr, seg.seg);
    }
    return arr;
}

void algorithm::run(std::string input_path) {

    // Load instance using the library
    auto data = SHPLoader::ReadShapeFileToPoint2D(input_path);
    std::cout << "Successfully loaded instance from: " << input_path << std::endl;
    std::cout << "Points: " << data.first.size() << std::endl;
    std::cout << "Polygons: " << data.second.size() << std::endl;
    std::vector<Segment_w_info> segments_w_info = build_segment_vector(data.first, data.second);
    std::vector<Segment> segments= segs_wo_info(segments_w_info);
    /*
    std::vector<Segment> segments;
    segments.emplace_back(Point(1.0,0.0), Point(1.0, 4.0));
    segments.emplace_back(Point(4.0, 0.0), Point(4.0, 4.0));
    segments.emplace_back(Point(2.0, 1.0), Point(3.0, 1.0));
    segments.emplace_back(Point(2.0, 3.0), Point(3.0, 3.0));
    segments.emplace_back(Point(2.0, 0.0), Point(3.0, 0.0));
    segments.emplace_back(Point(2.0, 4.0), Point(3.0, 4.0));
    segments.emplace_back(Point(3.0, 1.0), Point(3.0, 3.0));
    std::vector<Segment_w_info> segments_w_info;
    int id=0;
    for (auto seg : segments) {
        segments_w_info.emplace_back(Segment_w_info(seg, true, id));
        ++id;
    }*/

    segments_to_svg(segments,"original.svg");
    auto extended_segs = ray_shoot_intersection(segments_w_info);
    std::cout<<"Extended Segments size: "<<extended_segs.size()<<std::endl;
    segments_to_svg(segs_wo_info(extended_segs),"extended.svg");

    Arrangement arr = build_arrangement(extended_segs);

    std::cout << "Arrangement has "
          << arr.number_of_vertices() << " vertices, "
          << arr.number_of_edges() << " edges, "
          << arr.number_of_faces() << " faces." << std::endl;
}