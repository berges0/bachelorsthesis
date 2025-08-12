//
// Created by samuel-berge on 8/11/25.
//
#include "core/edge_extension.hpp"
#include <algorithm>

namespace EDGE_EXTENSION {

    const std::vector<Segment_w_info> &edge_extension(std::vector<Segment_w_info>& segments, std::string version) {
        if (version == "0") {
            return STANDARD::extension(segments);
        }
        if (version == "1") {
            return LIMITED::extension(segments);
        }
        else {
            throw std::runtime_error("Unknown version: " + version);
        }
        static std::vector<Segment_w_info> empty;
        return empty; // This should never be reached, but added to avoid compiler warnings.
    }

    namespace STANDARD {
        const std::vector<Segment_w_info> &extension(std::vector<Segment_w_info>& segments) {

            std::vector<Segment> just_segments = filter_segments(segments);
            Tree tree(just_segments.begin(), just_segments.end());
            tree.build();

            static std::vector<Segment_w_info> all_segments = segments;

            int count = 0;
            for (const Segment& segment : just_segments) {
                Line line(segment);
                Point p1 = segment.source();
                Point p2 = segment.target();
                //std::cout << p1.x() << ", "<< p1.y() << " / " << p2.x() <<", "<< p2.y() << "\n";

                Vector forward = p1 - p2;
                Vector backward = p2 - p1;

                if (segments[count].shoot_source) {
                    auto hit_forward = first_intersection(p1, forward, tree, &segment);
                    if (hit_forward) {
                        all_segments.emplace_back(Segment_w_info(Segment(p1, *hit_forward), false, -1, false,false));
                        //std::cout<< "EMPLACED " << all_segments.back().seg.source().x() << " " << all_segments.back().seg.source().y()<<", "<< all_segments.back().seg.target().x()<<" "<<all_segments.back().seg.target().y()<<std::endl;
                    }

                }
                if (segments[count].shoot_target) {
                    auto hit_backward = first_intersection(p2, backward, tree, &segment);
                    if (hit_backward) {
                        all_segments.emplace_back(Segment_w_info(Segment(p2, *hit_backward),false, -1,false, false));
                        //std::cout<< "EMPLACED " << all_segments.back().seg.source().x() << " " << all_segments.back().seg.source().y()<<", "<< all_segments.back().seg.target().x()<<" "<<all_segments.back().seg.target().y()<<std::endl;
                    }
                }
                count++;
            }
        return all_segments;
        }

    }

    namespace LIMITED{
        const std::vector<Segment_w_info> &extension(std::vector<Segment_w_info>& segments) {
            return segments;
            post_process(segments);
        }

        void post_process(std::vector<Segment_w_info>& segments) {
        }
    }

    std::optional<Point> first_intersection(const Point& origin,
    const Vector& direction,
    const Tree& tree,
    const Segment* self_segment) {

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

    std::vector<Segment> filter_segments(const std::vector<Segment_w_info>& segments) {

        std::vector<Segment> just_segments;
        for (const auto& seg : segments) {
            just_segments.emplace_back(seg.seg);
        }
        return just_segments;
    }


    void add_outer_box(std::vector<Segment_w_info>& segments) {

        double minX = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        double maxY = std::numeric_limits<double>::lowest();

        for (const auto& seg : segments) {
            minX = std::min(minX,std::min(CGAL::to_double(seg.seg.source().x()), CGAL::to_double(seg.seg.target().x())));
            maxX = std::max(maxX, std::max(CGAL::to_double(seg.seg.source().x()), CGAL::to_double(seg.seg.target().x())));
            minY = std::min(minY, std::min(CGAL::to_double(seg.seg.source().y()), CGAL::to_double(seg.seg.target().y())));
            maxY = std::max(maxY, std::max(CGAL::to_double(seg.seg.source().y()), CGAL::to_double(seg.seg.target().y())));
        }
        // Create the outer box segments
        int id = segments.size();
        double set_off_x = 0.01*(maxX - minX);
        double set_off_y = 0.01*(maxY - minY);
        Point blc(minX-set_off_x, minY-set_off_y);
        Point tlc(minX-set_off_x, maxY+set_off_y);
        Point trc(maxX+set_off_x, maxY+set_off_y);
        Point brc(maxX+set_off_x, maxY-set_off_y);

        Segment_w_info edge1(Segment(blc, tlc), false, id++,false,false);
        Segment_w_info edge2(Segment(tlc, trc), false, id++,false,false);
        Segment_w_info edge3(Segment(trc, brc), false, id++,false,false);
        Segment_w_info edge4(Segment(brc, blc), false, id++,false,false);
        segments.push_back(edge1);
        segments.push_back(edge2);
        segments.push_back(edge3);
        segments.push_back(edge4);
    }

}
