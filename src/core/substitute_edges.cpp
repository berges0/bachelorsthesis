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
        auto relinked = order_endpoints_by_orthogonal_projection(EDGE_EXTENSION::filter_segments(group));
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
        }
    }
    IO_FUNCTIONS::SVG::segments_to_svg(test, "testttttt.svg");
}


static Kernel::FT dot(const Vector& a, const Vector& b) {
    return a.x()*b.x() + a.y()*b.y();
}
struct EndpointParam {
    Point p;
    Kernel::FT    s;   // scalar parameter along the averaged line (c + s * dir)
};

std::vector<Point> order_endpoints_by_orthogonal_projection(const std::vector<Segment>& segs)
{
    // Gather endpoints and build a predominant direction by sign-aligned summation.
    std::vector<Point> endpoints;
    endpoints.reserve(segs.size() * 2);

    Vector sum(Kernel::FT(0), Kernel::FT(0));
    Vector ref;               // reference to fix sign
    bool   have_ref = false;

    for (const auto& s : segs) {
        endpoints.push_back(s.source());
        endpoints.push_back(s.target());

        Vector v = s.to_vector();
        if (v == CGAL::NULL_VECTOR) continue;

        if (!have_ref) { ref = v; have_ref = true; }
        // Flip direction if mostly opposite to the reference.
        if (dot(v, ref) < Kernel::FT(0)) v = -v;
        sum = sum + v;
    }
    if (!have_ref || sum == CGAL::NULL_VECTOR) {
        // Fallback: nothing to order by → return in input order.
        return endpoints;
    }

    // Define the averaged line: through the centroid c with direction "sum".
    //   (No normalization needed; we’ll divide by |sum|^2 when projecting.)
    Kernel::FT cx(0), cy(0);
    for (const auto& p : endpoints) { cx += p.x(); cy += p.y(); }
    const Kernel::FT n = Kernel::FT( (int)endpoints.size() );
    const Point c( cx / n, cy / n );

    const Kernel::FT denom = dot(sum, sum); // = |sum|^2  (positive)

    // Orthogonally project each endpoint onto that line and record its parameter s.
    std::vector<EndpointParam> proj;
    proj.reserve(endpoints.size());
    for (const auto& p : endpoints) {
        Vector cp = p - c;                 // vector from c to p
        Kernel::FT s = dot(cp, sum) / denom;       // parameter along the line: c + s * sum
        proj.push_back({p, s});
    }

    // Sort by the parameter along the line (left → right along the main direction).
    std::sort(proj.begin(), proj.end(),
              [](const EndpointParam& a, const EndpointParam& b){ return a.s < b.s; });

    // Emit ordered endpoints.
    std::vector<Point> ordered;
    ordered.reserve(proj.size());
    for (const auto& e : proj) ordered.push_back(e.p);
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