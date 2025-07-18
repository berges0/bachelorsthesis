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
        if(n<2) continue;

        if(n==2){
            segments.emplace_back(Segment_w_info(Segment( Point(points[poly[0]].x, points[poly[0]].y), Point(points[poly[1]].x, points[poly[1]].y)), true, id++, true, true));
        }

        bool last_shoot;
        if (CGAL::orientation(Point(points[poly[n-1]].x, points[poly[n-1]].y), Point(points[poly[0]].x, points[poly[0]].y), Point(points[poly[1]].x, points[poly[1]].y)) == CGAL::RIGHT_TURN) {
            last_shoot=false;

        }
        else if (CGAL::orientation(Point(points[poly[n-1]].x, points[poly[n-1]].y), Point(points[poly[0]].x, points[poly[0]].y), Point(points[poly[1]].x, points[poly[1]].y)) == CGAL::LEFT_TURN) {
            last_shoot=true;
        }
        else{
            throw std::runtime_error("instance error: straight line consists of more than one edge");
        }

        bool forward_shoot;

        for(int i=0; i<n; i++){

            Point pt1(points[poly[i]].x, points[poly[i]].y);
            Point pt2(points[poly[((i+1) % n)]].x, points[poly[((i+1) % n)]].y);
            Point pt3(points[poly[((i+2) % n)]].x, points[poly[((i+2) % n)]].y);

            //UPDATE LIMITS
            if (points[poly[i]].x < minX) minX = points[poly[i]].x;
            if (points[poly[i]].x > maxX) maxX = points[poly[i]].x;
            if (points[poly[i]].y < minY) minY = points[poly[i]].y;
            if (points[poly[i]].y > maxY) maxY = points[poly[i]].y;

            if (CGAL::orientation(pt1, pt2, pt3) == CGAL::RIGHT_TURN) {
                forward_shoot=false;

            } else if (CGAL::orientation(pt1, pt2, pt3) == CGAL::LEFT_TURN) {
                forward_shoot=true;
            } else{std::cerr<<"ERROR";}

            segments.emplace_back(Segment_w_info(Segment(pt1,pt2), true, id++, last_shoot, forward_shoot));
            last_shoot=forward_shoot;
        }
    }
    std::pair<Point, Point> limits = {Point(minX, minY), Point(maxX, maxY)};
    return {segments, limits};
}

std::vector<Segment> algorithm::segs_wo_info(const std::vector<Segment_w_info>& segments) {
    std::vector<Segment> just_segments;
    for (const auto& seg : segments) {
        just_segments.emplace_back(seg.seg);
    }
    return just_segments;
}
