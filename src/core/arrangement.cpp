//
// Created by samuel-berge on 8/11/25.
//
#include "core/arrangement.hpp"

namespace ARRANGEMENT {

    Arrangement build_arrangement(const std::vector<Segment_w_info>& segments) {
        Arrangement arr;
        std::vector<Curve> curves;
        int next_id = 0;
        int count11 = 0;

        for (const auto& seg : segments) {
            Curve c(seg.seg, seg.from_poly);
            //std::cout << "Added segment from" << seg.from_poly << " polygon " << std::endl;
            curves.push_back(c);
        }

        CGAL::insert(arr, curves.begin(), curves.end());
        //std::cout<< "Arrangement has "<< arr.number_of_edges() << std::endl;
        add_edge_data(arr);
        add_face_data(arr);

        return arr;
    }

    void add_edge_data(Arrangement &arr) {

        for(auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
            auto& edge = *eit;
            auto curve = edge.curve();
            const auto& unique_list = curve.data();
            bool from_polygon = (*unique_list.begin())->data();
            edge.set_data(HalfedgeData{from_polygon});
            edge.twin()->set_data(HalfedgeData{from_polygon});
        }
    }

    void add_face_data(Arrangement &arr) {

        size_t face_count = 0;
        Polygon_2 poly;

        for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
            poly.clear();
            bool pure_poly;
            if (fit->is_unbounded()) {
                fit->set_data(FaceData(-1, false, DBL_MAX));
                continue;
            }
            pure_poly=true;
            // outer boundary exists
            auto circ = fit->outer_ccb();
            auto curr = circ;
            do {
                const auto& source = curr->source()->point();
                const auto& target = curr->target()->point();
                poly.push_back(source);
                if (!(curr->data().frompoly)) {
                    pure_poly = false;
                }
                ++curr;
            } while (curr != circ);

            assert(poly.is_simple());

            double area = CGAL::to_double(poly.area());
            fit->set_data(FaceData(face_count, pure_poly, area));
            ++face_count;
        }
    }


    Arrangement build_arrangement_relinked(const std::vector<Segment_w_info>& segments, RTree &rtree, const std::vector<PWH> &polygonswh) {
        Arrangement arr;
        std::vector<Curve> curves;
        int next_id = 0;
        int count11 = 0;

        for (const auto& seg : segments) {
            Curve c(seg.seg, seg.from_poly);
            //std::cout << "Added segment from" << seg.from_poly << " polygon " << std::endl;
            curves.push_back(c);
        }

        CGAL::insert(arr, curves.begin(), curves.end());
        //std::cout<< "Arrangement has "<< arr.number_of_edges() << std::endl;
        add_edge_data(arr);
        add_face_data_relinked(arr, rtree, polygonswh);

        return arr;
    }

    void add_face_data_relinked(Arrangement &arr, RTree &rtree, const std::vector<PWH> &polygonswh) {
        size_t face_count = 0;
        Polygon_2 poly;

        for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
            poly.clear();
            bool pure_poly=false;
            bool could_be_poly=false;
            if (fit->is_unbounded()) {
                fit->set_data(FaceData(-1, false, DBL_MAX));
                continue;
            }
            // outer boundary exists
            auto circ = fit->outer_ccb();
            auto curr = circ;
            do {
                const auto& source = curr->source()->point();
                const auto& target = curr->target()->point();
                poly.push_back(source);
                if ((curr->data().frompoly)) {
                    could_be_poly = true;
                }
                ++curr;
            } while (curr != circ);

            if (could_be_poly) {
                auto a = curr->source()->point();
                auto b = curr->target()->point();

                // midpoint
                Point m = CGAL::midpoint(a,b);

                // direction vector
                Vector d = b - a;

                // left normal (-dy, dx)
                Vector n(-d.y(), d.x());

                // normalize
                double len = std::sqrt(CGAL::to_double(n.squared_length()));
                n = n / len;   // unit vector to the left

                // epsilon shift
                double eps = 0.00001; // pick your offset size
                Point shifted = m + eps * n;
                pure_poly = test_in_poly(shifted, rtree, polygonswh);
            }

            assert(poly.is_simple());

            double area = CGAL::to_double(poly.area());
            fit->set_data(FaceData(face_count, pure_poly, area));
            ++face_count;
        }
    }

bool test_in_poly(const Point &qp, const RTree &rtree, const std::vector<PWH> &polys) {
        BBox qbox(BPoint(CGAL::to_double(qp.x()), CGAL::to_double(qp.y())),
        BPoint(CGAL::to_double(qp.x()), CGAL::to_double(qp.y())));

        // --- Candidate lookup + exact test ---
        std::vector<RItem> hits;
        rtree.query(bgi::intersects(qbox), std::back_inserter(hits));

        bool inside_any = false;
        for (const auto& it : hits) {
            const Polygon_2& outer = polys[it.second].outer_boundary(); // <- use outer polygon
            auto side = outer.bounded_side(qp);
            if (side != CGAL::ON_UNBOUNDED_SIDE) { // counts boundary as inside
                return true;
            }
        }
        return false;
    }
}
