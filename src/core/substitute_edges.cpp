//
// Created by samuel-berge on 8/30/25.
//

#include "core/substitute_edges.hpp"

#include "core/edge_extension.hpp"
#include "io/io_functions.hpp"

namespace SUBSTITUTE_EDGES {
void relink_edges(std::vector<std::vector<Segment_w_info>> &segments, Tree &tree) {
    auto *skip = &segments[0];
    std::vector <Segment> test;
    for (auto &group : segments) {
        if ((&group == skip)) continue;
        if (group.size() == 1){std::cout << "SOMETHING BAD HAPPENED"<< std::endl;}
        //IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(group), std::to_string(k++)+"_group_b.svg");
        auto relinked = order_endpoints_along_main_dir(EDGE_EXTENSION::filter_segments(group));
        int count = 0;
        for (int i = 0; i < relinked.size() - 1; ++i) {
            Segment seg = Segment(relinked[i], relinked[i + 1]);
            if (seg.squared_length() == 0)continue;

            Point mid_point = CGAL::midpoint(seg.source(), seg.target());
            if (tree.any_intersection(mid_point)) {
                test.push_back(seg);
                continue;
            };
            segments[0].emplace_back(seg, false, -1,
                -1, false, false, true);
            count++;
        }
    }
    IO_FUNCTIONS::SVG::segments_to_svg(test, "testttttt.svg");
    throw std::runtime_error("STOP");
}

struct EndpointAlongDir {
    Point p;
    double t; // projection parameter
};

std::vector<Point> order_endpoints_along_main_dir(const std::vector<Segment>& segs) {
    // 1) Gather endpoints and estimate the predominant direction.
    std::vector<Point> endpoints;
    endpoints.reserve(segs.size() * 2);

    double sx = 0.0, sy = 0.0;      // accumulated (signed) unit directions
    double refx = 0.0, refy = 0.0;  // reference direction to fix sign flips
    bool   have_ref = false;

    for (const auto& s : segs) {
        endpoints.push_back(s.source());
        endpoints.push_back(s.target());

        Vector v = s.to_vector();
        double vx = CGAL::to_double(v.x());
        double vy = CGAL::to_double(v.y());
        double L  = std::hypot(vx, vy);
        if (L == 0.0) continue; // skip degenerate segments

        double ux = vx / L, uy = vy / L; // unit direction of this segment

        if (!have_ref) { refx = ux; refy = uy; have_ref = true; }
        // Flip if pointing mostly the opposite way, to keep a consistent sign.
        if (ux*refx + uy*refy < 0.0) { ux = -ux; uy = -uy; }

        sx += ux; sy += uy;
    }

    // Fallback if everything was degenerate
    if (!have_ref) return endpoints;

    // Normalize the accumulated direction
    double sL = std::hypot(sx, sy);
    if (sL == 0.0) { sx = refx; sy = refy; sL = std::hypot(sx, sy); }
    double dx = sx / sL, dy = sy / sL; // predominant unit direction

    // 2) Choose an origin (centroid of endpoints works well for numerical stability).
    double cx = 0.0, cy = 0.0;
    for (const auto& p : endpoints) { cx += CGAL::to_double(p.x()); cy += CGAL::to_double(p.y()); }
    if (!endpoints.empty()) { cx /= endpoints.size(); cy /= endpoints.size(); }

    // 3) Project endpoints onto the direction and sort by the scalar parameter t.
    std::vector<EndpointAlongDir> proj;
    proj.reserve(endpoints.size());
    for (const auto& p : endpoints) {
        double px = CGAL::to_double(p.x());
        double py = CGAL::to_double(p.y());
        double t  = (px - cx)*dx + (py - cy)*dy; // signed distance along the direction
        proj.push_back({p, t});
    }
    std::sort(proj.begin(), proj.end(),
              [](const EndpointAlongDir& a, const EndpointAlongDir& b){ return a.t < b.t; });

    // 4) Return endpoints in order (left→right along predominant direction).
    std::vector<Point> ordered;
    ordered.reserve(proj.size());
    for (auto& e : proj) ordered.push_back(e.p);
    return ordered;
}

RTree build_r_tree(const std::vector<PWH> &polys) {
    std::vector<RItem> items;
    items.reserve((int)polys.size());
    for (int i = 0; i < (int)polys.size(); ++i) {
        CGAL::Bbox_2 b = polys[i].bbox();
        items.emplace_back(BBox(BPoint(b.xmin(), b.ymin()), BPoint(b.xmax(), b.ymax())), i);
    }
    RTree rtree(items.begin(), items.end());
    return rtree;
}

void connect_outer_points(std::vector<std::vector<Segment_w_info>> &segments) {
    auto *skip = &segments[0];
    //int k = 0;
    for (auto &group : segments) {
        if ((&group == skip)||group.size() <= 1) continue;
        auto order = order_endpoints_along_main_dir(EDGE_EXTENSION::filter_segments(group));
        std::vector<Segment_w_info> segs;
        //IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(group), "group_" + std::to_string(k) + ".svg");
        //std::cout<<"Distance "<<std::sqrt(CGAL::to_double(CGAL::squared_distance(group[0].seg, group[group.size()-1].seg)))<<std::endl;
        Segment seg = Segment(order[0], order[order.size()-1]);
        auto control = EDGE_EXTENSION::filter_segments(group);
        control.push_back(seg);
        //IO_FUNCTIONS::SVG::segments_to_svg(control, "group_" + std::to_string(k++) + "with_connection.svg");
        if (seg.squared_length() == 0) continue;
        segments[0].emplace_back(seg, false, -1,
            -1, false, false, true);
    }
}

void post_prune(std::vector<Segment_w_info> &segments) {
    std::vector<Segment> just_segments = EDGE_EXTENSION::filter_segments(segments);
    Tree tree(just_segments.begin(), just_segments.end());
    tree.build();
    tree.accelerate_distance_queries();


    int seg_id = segments.size();

    std::vector<Segment_w_info> pruned_segments;
    for (size_t i = 0; i < segments.size(); ++i) {
        if (!segments[i].from_poly) {
            Point P1 = segments[i].seg.source();
            Point P2 = segments[i].seg.target();


            Vector dir_forward = P1 - P2;
            Vector dir_backward = P2 - P1;

            std::optional<Point> hit_forward;
            std::optional<Point> hit_backward;

            if (tree.number_of_intersected_primitives(P1)>1) {
                hit_forward = P1;
            }
            else {
                hit_forward = EDGE_EXTENSION::first_intersection(P1, dir_forward, tree, &just_segments[i]);
            }
            if (tree.number_of_intersected_primitives(P2)>1) {
                hit_backward = P2;
            }
            else {
                hit_backward = EDGE_EXTENSION::first_intersection(P2, dir_backward, tree, &just_segments[i]);
            }
            //two valid intersections, not the same point and just a shorting of original edge
            if (hit_backward && hit_forward && EDGE_EXTENSION::get_distance(P2,*hit_backward)<=EDGE_EXTENSION::get_distance(P1,P2)&&
                EDGE_EXTENSION::get_distance(P1,*hit_forward)<=EDGE_EXTENSION::get_distance(P1,P2) && *hit_forward != *hit_backward){
                pruned_segments.emplace_back(Segment_w_info(Segment(*hit_backward, *hit_forward), segments[i].from_poly, seg_id++,
                    segments[i].poly_id, segments[i].shoot_source,segments[i].shoot_target, segments[i].replacing_edge));
                }
        }
        else {
            pruned_segments.emplace_back(Segment_w_info(segments[i].seg, segments[i].from_poly, seg_id++, segments[i].poly_id,
                segments[i].shoot_source,segments[i].shoot_target, segments[i].replacing_edge));
        }
    }
    segments = pruned_segments;
}

}