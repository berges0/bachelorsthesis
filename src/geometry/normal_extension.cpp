#include "algorithm.hpp"

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

    static std::vector<Segment_w_info> all_segments = segments;


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
