#include "algorithm.hpp"

Arrangement &algorithm::build_arrangement(std::vector<Segment_w_info> segments) {
    static Arrangement arr;
    std::vector<Curve> curves;
    int next_id = 0;
    for (const auto& seg : segments) {
        Curve c(seg.seg, seg.from_poly);
        //std::cout << "Added segment from" << seg.from_poly << " polygon " << std::endl;
        curves.push_back(c);
    }
    CGAL::insert(arr, curves.begin(), curves.end());
    for(auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
        auto& edge = *eit;
        auto curve = edge.curve();
        const auto& unique_list = curve.data();
        bool from_polygon = (*unique_list.begin())->data(); 
        edge.set_data(HalfedgeData{from_polygon});
        edge.twin()->set_data(HalfedgeData{from_polygon});
    }
    size_t face_count = 0;
    Polygon_2 poly;
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        poly.clear();
        bool pure_poly;
        if (fit->is_unbounded()) {
            fit->set_data(FaceData(0, pure_poly, DBL_MAX));
            continue;
        }
        else {
            pure_poly=true;
            // outer boundary exists
            auto circ = fit->outer_ccb();
            auto curr = circ;
            do {
                const auto& source = curr->source()->point();
                const auto& target = curr->target()->point();
                poly.push_back(source);
                //std::cout << "  Edge: (" << source << ") -> (" << target << ")\n";
                if (!(*curr).data().frompoly) {
                    pure_poly = false;
                }
                ++curr;
            } while (curr != circ);
        }
        if (!poly.is_simple()) {
    //std::cerr << "Polygon for face " << face_count << " is NOT simple!" << std::endl;
}
        double area = CGAL::to_double(poly.area());
        //std::cout<<"FACE "<<fit->data().id<<" area: "<<area<<std::endl;
        fit->set_data(FaceData(face_count, pure_poly, area));
        // Access data field here if your faces have one:
        // std::cout << "Data: " << face.data() << "\n";
        ++face_count;
    }
    return arr;
}

