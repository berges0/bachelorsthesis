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

}
