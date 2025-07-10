#include "algorithm.hpp"

std::vector<Segment> algorithm::segs_wo_info(const std::vector<Segment_w_info>& segments) {
    std::vector<Segment> just_segments;
    for (const auto& seg : segments) {
        just_segments.emplace_back(seg.seg);
    }
    return just_segments;
}


std::pair<std::vector<Segment_w_info>, std::pair<Point, Point>> algorithm::build_segment_vector(const std::vector<SHPLoader::Point2D>& points,
                    const std::vector<std::vector<int>>& polygons) {
    std::vector<Segment_w_info> segments;
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

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

    return {segments, limits};
}

void algorithm::add_box(std::vector<Segment_w_info> &segs, std::pair<Point, Point> limits) {
    int id = segs.size();
    Point blc(limits.first.x(), limits.first.y());
    Point tlc(limits.first.x(), limits.second.y());
    Point trc(limits.second.x(), limits.second.y());
    Point brc(limits.second.x(), limits.first.y());

    Segment_w_info edge1(Segment(blc, tlc), false, id++);
    Segment_w_info edge2(Segment(tlc, trc), false, id++);
    Segment_w_info edge3(Segment(trc, brc), false, id++);
    Segment_w_info edge4(Segment(brc, blc), false, id++);
    segs.push_back(edge1);
    segs.push_back(edge2);
    segs.push_back(edge3);
    segs.push_back(edge4);
}



std::optional<Point> algorithm::nearest_intersection_in_direction(const Point& origin,
    const Vector& direction,
    const Tree& tree,
    const Segment* self_segment){
    Ray ray(origin, origin + direction);
    SkipSegment skip(self_segment);

    auto result = tree.first_intersection(ray,skip);
    if (result) {
        //std::cout << "Intersection found\n";
        // unpack intersection as before, return point
    } else {
        //std::cout << "No intersection found from origin: "
                  //<< CGAL::to_double(origin.x()) << ", "
                  //<< CGAL::to_double(origin.y()) << "\n";
    }
    if (result) {
        if (const Point* ipoint = std::get_if<Point>(&(result->first))) {
            if (*ipoint != self_segment->source() && *ipoint != self_segment->target()) {
                //std::cout<<"and it is not itself but "<<origin<< ", "<< *ipoint <<std::endl;
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
        //std::cout << p1.x() << ", "<< p1.y() << " / " << p2.x() <<", "<< p2.y() << "\n";

        Vector forward = p1 - p2;
        Vector backward = p2 - p1;

        auto hit_forward = nearest_intersection_in_direction(p1, forward, tree, &segment);
        auto hit_backward = nearest_intersection_in_direction(p2, backward, tree, &segment);

        if (hit_forward) {
            all_segments.emplace_back(Segment_w_info(Segment(p1, *hit_forward), false, -1));
            //std::cout<< "EMPLACED " << all_segments.back().seg.source().x() << " " << all_segments.back().seg.source().y()<<", "<< all_segments.back().seg.target().x()<<" "<<all_segments.back().seg.target().y()<<std::endl;
        }
        if (hit_backward) {
            all_segments.emplace_back(Segment_w_info(Segment(p2, *hit_backward),false, -1));
            //std::cout<< "EMPLACED " << all_segments.back().seg.source().x() << " " << all_segments.back().seg.source().y()<<", "<< all_segments.back().seg.target().x()<<" "<<all_segments.back().seg.target().y()<<std::endl;
        }
    }
    return all_segments;
}


double algorithm::calculate_length(const Arrangement::Halfedge_const_handle& he) {
    auto p1 = he->source()->point();
    auto p2 = he->target()->point();

    auto dx = p2.x() - p1.x();
    auto dy = p2.y() - p1.y();

    auto squared_length = dx*dx + dy*dy;
    return std::sqrt(CGAL::to_double(squared_length));
}