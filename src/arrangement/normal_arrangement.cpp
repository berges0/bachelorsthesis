#include "algorithm.hpp"

Arrangement &algorithm::build_arrangement(std::vector<Segment_w_info> segments) {
    static Arrangement arr;
    std::vector<Curve> curves;
    int next_id = 0;
    int count11 = 0;

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
        if (from_polygon) {
            count11++;
        }
        edge.set_data(HalfedgeData{from_polygon});
        edge.twin()->set_data(HalfedgeData{from_polygon});
    }

    size_t face_count = 0;
    Polygon_2 poly;
    int countnew=0;
    int count5=0;

    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        for (auto edge = fit->outer_ccb(); edge != fit->outer_ccb(); ++edge) {
            if (edge->data().frompoly) {
                fit->set_data(FaceData(face_count, true));
                break; // only need to set once
            }
        }
        poly.clear();
        bool pure_poly;
        if (fit->is_unbounded()) {
            fit->set_data(FaceData(-1, false, DBL_MAX));
            count5++;
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

            //std::cout << "  Edge: (" << source << ") -> (" << target << ")\n";
            if (!(curr->data().frompoly)) {

                pure_poly = false;
            }
            ++curr;
        } while (curr != circ);

            countnew++;

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
    int count10 = 0;
    std::cout <<"There are " << countnew << " polyfaces " <<std::endl;
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (fit->is_unbounded()) {
            continue; // skip unbounded faces
        }
        if (fit->data().belongs_to_poly) {
            count10++;
        }
    }
    std::vector<Segment> segs;
    for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
        if (eit->data().frompoly){
            assert(eit->twin()->data().frompoly);
            const auto& source = eit->source()->point();
            const auto& target = eit->target()->point();
            segs.emplace_back(Segment(source, target));
        }
    }
    segments_to_svg(segs, "checkpoint.svg");
    return arr;
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

