#include "utils/metric_simplification.hpp"

namespace METRIC_SIMPLIFICATION {

// ----------------- Hilfsstrukturen & Utilities -----------------------------

struct Hit {
    bool ok = false;
    Point p;
    const Segment* seg = nullptr; // zeigt in Segmentcontainer (Lebensdauer beachten)
    Kernel::FT sqdist = Kernel::FT(0);
};

static inline Vector left_normal(const Vector& v)  { return Vector(-v.y(),  v.x()); }
static inline Vector right_normal(const Vector& v) { return Vector( v.y(), -v.x()); }

static inline double seg_length_d(const Segment& s) {
    return std::sqrt(CGAL::to_double(s.squared_length()));
}

static inline double total_length_d(std::vector<Segment>& segs) {
    Kernel::FT L(0);
    for (const auto& s : segs) L += std::sqrt(CGAL::to_double(s.squared_length()));
    return CGAL::to_double(L);
}

// Winkelbetrag zweier Segmente in Grad (0..90)
static inline double angle_deg_abs(const Segment& a, const Segment& b) {
    Vector va = a.to_vector(), vb = b.to_vector();
    Kernel::FT la2 = va.squared_length();
    Kernel::FT lb2 = vb.squared_length();
    if (la2 == Kernel::FT(0) || lb2 == Kernel::FT(0)) return 0.0;

    Kernel::FT dot = va.x()*vb.x() + va.y()*vb.y();
    double num = std::fabs(CGAL::to_double(dot));
    double den = std::sqrt(CGAL::to_double(la2) * CGAL::to_double(lb2));
    double c = (den > 0.0) ? std::min(1.0, num/den) : 1.0;
    double rad = std::acos(c);
    return rad * 180.0 / M_PI;
}

// Endpunkt eines begrenzten Normalsegments: origin + dir * (max_dist / |dir|)
static inline bool make_bounded_endpoint(const Point& origin,
                                         const Vector& dir,
                                         double max_dist,
                                         Point& out_end)
{
    Kernel::FT len2 = dir.squared_length();
    if (len2 == Kernel::FT(0)) return false;

    double dir_len = std::sqrt(CGAL::to_double(len2));
    if (dir_len == 0.0) return false;

    double scale = max_dist / dir_len;
    double dx = CGAL::to_double(dir.x()) * scale;
    double dy = CGAL::to_double(dir.y()) * scale;

    out_end = Point(origin.x() + Kernel::FT(dx),
                    origin.y() + Kernel::FT(dy));
    return true;
}

// 2-Arg Overload von AABB_tree::all_intersections
static Hit first_hit_bounded(const Point& origin,
                             const Vector& dir,
                             double max_dist,
                             const Tree& tree)
{
    Hit h{};
    Point end;
    if (!make_bounded_endpoint(origin, dir, max_dist, end)) return h;

    Segment query(origin, end);

    std::vector<typename Tree::Object_and_primitive_id> hits;
    tree.all_intersections(query, std::back_inserter(hits));

    if (hits.empty()) return h;

    Kernel::FT best = Kernel::FT(-1);
    const Segment* best_seg = nullptr;
    Point best_p;

    for (const auto& op : hits) {
        Point ip;
        Segment iseg;
        if (CGAL::assign(ip, op.first)) {
            Kernel::FT d = CGAL::squared_distance(origin, ip);
            if (best < Kernel::FT(0) || d < best) {
                best = d; best_p = ip; best_seg = &*op.second;
            }
        } else if (CGAL::assign(iseg, op.first)) {
            // Kollinearüberlappung: nimm näheren Endpunkt
            Kernel::FT d1 = CGAL::squared_distance(origin, iseg.source());
            Kernel::FT d2 = CGAL::squared_distance(origin, iseg.target());
            Point ip2 = (d1 <= d2) ? iseg.source() : iseg.target();
            Kernel::FT d  = (d1 <= d2) ? d1 : d2;
            if (best < Kernel::FT(0) || d < best) {
                best = d; best_p = ip2; best_seg = &*op.second;
            }
        }
    }

    if (best_seg) {
        h.ok = true; h.p = best_p; h.seg = best_seg; h.sqdist = best;
    }
    return h;
}

static inline bool unit_vector(const Vector& v, Vector& u_out) {
    Kernel::FT l2 = v.squared_length();
    if (l2 == Kernel::FT(0)) return false;
    double invl = 1.0 / std::sqrt(CGAL::to_double(l2));
    u_out = Vector(v.x() * Kernel::FT(invl), v.y() * Kernel::FT(invl));
    return true;
}

static inline Point advance_point(const Point& p, const Vector& dir_unit, double dist) {
    return Point(p.x() + Kernel::FT(CGAL::to_double(dir_unit.x()) * dist),
                 p.y() + Kernel::FT(CGAL::to_double(dir_unit.y()) * dist));
}

// ----------------- Segment-Extraktion -------------------------------------

static inline void append_segments_from_polygon(const Polygon_2& P, std::vector<Segment>& out) {
    if (P.is_empty()) return;
    out.reserve(out.size() + P.size());
    for (auto it = P.edges_begin(); it != P.edges_end(); ++it) out.push_back(*it);
}

static inline void append_segments_from_pwh(const PWH& pwh, std::vector<Segment>& out) {
    append_segments_from_polygon(pwh.outer_boundary(), out);
    for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
        append_segments_from_polygon(*hit, out);
    }
}

static inline std::vector<Segment> segments_from(const std::vector<PWH>& polys) {
    std::vector<Segment> segs;
    for (const auto& P : polys) append_segments_from_pwh(P, segs);
    return segs;
}

// gemeinsame Bounding-Box-Diagonale für A und B
static inline double bbox_diag(const std::vector<Segment>& A,
                               const std::vector<Segment>& B)
{
    double minx= std::numeric_limits<double>::infinity(),
           miny= std::numeric_limits<double>::infinity(),
           maxx=-std::numeric_limits<double>::infinity(),
           maxy=-std::numeric_limits<double>::infinity();
    auto upd=[&](const Segment& s){
        auto u=s.source(), v=s.target();
        double ux=CGAL::to_double(u.x()), uy=CGAL::to_double(u.y());
        double vx=CGAL::to_double(v.x()), vy=CGAL::to_double(v.y());
        minx=std::min(minx,std::min(ux,vx));
        miny=std::min(miny,std::min(uy,vy));
        maxx=std::max(maxx,std::max(ux,vx));
        maxy=std::max(maxy,std::max(uy,vy));
    };
    for (const auto& s:A) upd(s);
    for (const auto& s:B) upd(s);
    double dx=maxx-minx, dy=maxy-miny;
    return std::sqrt(dx*dx+dy*dy);
}

// ----------------- gerichtete Kernroutine: Ref -> Cmp (ohne/mit Viz) ------

static inline double directed_deviation_segments_core(std::vector<Segment>& Ref,
                                                      std::vector<Segment>& Cmp,
                                                      bool penalize_miss,
                                                      const RangeParams& rp,
                                                      const VizParams* vz,               // optional
                                                      std::vector<Segment>* viz_out)     // optional
{
    Tree tree(Cmp.begin(), Cmp.end()); // wichtig: iterator (nicht const_iterator)
    tree.accelerate_distance_queries();

    const double Ltot = total_length_d(Ref);
    if (Ltot == 0.0) return 0.0;

    // globale Skala aus Bounding Box
    const double diag = bbox_diag(Ref, Cmp);
    const double r_min_global = rp.min_factor * diag;
    const double r_max_global = rp.max_factor * diag;

    double acc = 0.0;

    for (const auto& e : Ref) {
        const double len = seg_length_d(e);
        if (len == 0.0) continue;

        // Reichweite r nach Modus wählen
        double r = 0.0;
        if (rp.use_fixed_range) {
            r = (rp.fixed_range_factor > 0.0) ? (rp.fixed_range_factor * diag)
                                              :  rp.fixed_range_abs;
            // optional klemmen
            if (rp.max_factor > 0.0) r = std::min(r, r_max_global);
            if (rp.min_factor > 0.0) r = std::max(r, r_min_global);
        } else {
            r = rp.alpha_len * len;               // abhängig von Kantenlänge
            r = std::max(r, r_min_global);
            r = std::min(r, r_max_global);
        }

        const Point mid(
            (e.source().x() + e.target().x())/Kernel::FT(2),
            (e.source().y() + e.target().y())/Kernel::FT(2)
        );

        const Vector v  = e.to_vector();
        const Vector nL = left_normal(v);
        const Vector nR = right_normal(v);

        const auto hL = first_hit_bounded(mid, nL, r, tree);
        const auto hR = first_hit_bounded(mid, nR, r, tree);

        const Hit* best = nullptr;
        const Vector* best_dir = nullptr;

        if (hL.ok && hR.ok) {
            bool left_is_better = (CGAL::compare(hL.sqdist, hR.sqdist) != CGAL::LARGER);
            best = left_is_better ? &hL : &hR;
            best_dir = left_is_better ? &nL : &nR;
        } else if (hL.ok) { best = &hL; best_dir = &nL; }
        else if (hR.ok)   { best = &hR; best_dir = &nR; }

        if (!best) {
            if (penalize_miss) acc += (len / Ltot) * 90.0;
            if (vz && viz_out && vz->draw_misses) {
                Vector u;
                if (unit_vector(nL, u)) {
                    Point end = advance_point(mid, u, r);
                    viz_out->push_back(Segment(mid, end));
                }
            }
            continue;
        }

        const double ang = angle_deg_abs(e, *best->seg);
        acc += (len / Ltot) * ang;

        // Visualisierung: prozentuales Overshoot mit Mindestwert, niemals länger als r
        if (vz && viz_out) {
            double hit_dist   = std::sqrt(CGAL::to_double(best->sqdist));
            double overshoot  = std::max(vz->overshoot_factor * hit_dist,
                                         vz->min_overshoot_factor * diag);
            double draw_len   = std::min(r, overshoot);
            Vector dir_u;
            if (unit_vector(*best_dir, dir_u)) {
                Point end = advance_point(mid, dir_u, hit_dist + draw_len);
                viz_out->push_back(Segment(mid, end));
            }
        }
    }
    return acc;
}

static inline double directed_deviation_segments(std::vector<Segment>& Ref,
                                                 std::vector<Segment>& Cmp,
                                                 bool penalize_miss,
                                                 const RangeParams& rp)
{
    return directed_deviation_segments_core(Ref, Cmp, penalize_miss, rp, nullptr, nullptr);
}

static inline double directed_deviation_segments_viz(std::vector<Segment>& Ref,
                                                     std::vector<Segment>& Cmp,
                                                     bool penalize_miss,
                                                     const RangeParams& rp,
                                                     const VizParams& vz,
                                                     std::vector<Segment>& viz_out)
{
    viz_out.clear();
    return directed_deviation_segments_core(Ref, Cmp, penalize_miss, rp, &vz, &viz_out);
}

// ----------------- öffentliche gerichtete/symmetrische APIs ----------------

double weighted_angle_deviation_bounded(const std::vector<PWH>& Aset,
                                        const std::vector<PWH>& Bset,
                                        bool penalize_miss,
                                        const RangeParams& rp)
{
    std::vector<Segment> segsA = segments_from(Aset);
    std::vector<Segment> segsB = segments_from(Bset);
    return directed_deviation_segments(segsA, segsB, penalize_miss, rp);
}

double weighted_angle_deviation_bounded_symmetric(const std::vector<PWH>& Aset,
                                                  const std::vector<PWH>& Bset,
                                                  bool penalize_miss,
                                                  const RangeParams& rp,
                                                  CombineMode mode)
{
    std::vector<Segment> segsA = segments_from(Aset);
    std::vector<Segment> segsB = segments_from(Bset);

    const double dAB = directed_deviation_segments(segsA, segsB, penalize_miss, rp);
    const double dBA = directed_deviation_segments(segsB, segsA, penalize_miss, rp);

    const double LA = total_length_d(segsA);
    const double LB = total_length_d(segsB);

    switch (mode) {
        case CombineMode::LengthWeightedMean: {
            const double denom = LA + LB;
            if (denom == 0.0) return 0.0;
            return (LA * dAB + LB * dBA) / denom;
        }
        case CombineMode::ArithmeticMean:
            return 0.5 * (dAB + dBA);
        case CombineMode::MaxOfDirected:
            return std::max(dAB, dBA);
        default:
            return 0.5 * (dAB + dBA);
    }
}

// ----------------- öffentliche APIs MIT Visualisierung ---------------------

double weighted_angle_deviation_bounded_with_viz(
    const std::vector<PWH>& Aset,
    const std::vector<PWH>& Bset,
    bool penalize_miss,
    const RangeParams& rp,
    const VizParams& vz,
    std::vector<Segment>& out_viz_AB)
{
    std::vector<Segment> segsA = segments_from(Aset);
    std::vector<Segment> segsB = segments_from(Bset);
    return directed_deviation_segments_viz(segsA, segsB, penalize_miss, rp, vz, out_viz_AB);
}

double weighted_angle_deviation_bounded_symmetric_with_viz(
    const std::vector<PWH>& Aset,
    const std::vector<PWH>& Bset,
    bool penalize_miss,
    const RangeParams& rp,
    const VizParams& vz,
    std::vector<Segment>& out_viz_AB,
    std::vector<Segment>& out_viz_BA,
    CombineMode mode)
{
    std::vector<Segment> segsA = segments_from(Aset);
    std::vector<Segment> segsB = segments_from(Bset);

    const double dAB = directed_deviation_segments_viz(segsA, segsB, penalize_miss, rp, vz, out_viz_AB);
    const double dBA = directed_deviation_segments_viz(segsB, segsA, penalize_miss, rp, vz, out_viz_BA);

    const double LA = total_length_d(segsA);
    const double LB = total_length_d(segsB);

    switch (mode) {
        case CombineMode::LengthWeightedMean: {
            const double denom = LA + LB;
            if (denom == 0.0) return 0.0;
            return (LA * dAB + LB * dBA) / denom;
        }
        case CombineMode::ArithmeticMean:
            return 0.5 * (dAB + dBA);
        case CombineMode::MaxOfDirected:
            return std::max(dAB, dBA);
        default:
            return 0.5 * (dAB + dBA);
    }
}

} // namespace METRIC_SIMPLIFICATION
