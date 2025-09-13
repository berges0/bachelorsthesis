//
// Created by samuel-berge on 8/11/25.
//
#include "io/io_functions.hpp"
//
// Created by samuel-berge on 8/10/25.
//

namespace IO_FUNCTIONS {
//SVG FILES

namespace SHP {

void read(const std::string filename, std::vector<Segment_w_info>& segments, Logger &logger) {

    auto data = SHPLoader::ReadShapeFileToPoint2D(filename);
    auto points = data.first;
    auto polygons = data.second;
    logger.add("Number of polygons in input", polygons.size());

    int seg_id=0;
    int poly_id=0;

    for (const auto &poly:polygons) {

        int n = poly.size();
        if(n<2) continue;

        if(n==2){
            segments.emplace_back(Segment_w_info(Segment( Point(points[poly[0]].x, points[poly[0]].y), Point(points[poly[1]].x, points[poly[1]].y)), true, seg_id++, poly_id++, true, true));
            continue;
        }

        bool backward_shoot;
        if (CGAL::orientation(Point(points[poly[n-1]].x, points[poly[n-1]].y), Point(points[poly[0]].x, points[poly[0]].y), Point(points[poly[1]].x, points[poly[1]].y)) == CGAL::RIGHT_TURN) {
            backward_shoot=false;

        }
        else if (CGAL::orientation(Point(points[poly[n-1]].x, points[poly[n-1]].y), Point(points[poly[0]].x, points[poly[0]].y), Point(points[poly[1]].x, points[poly[1]].y)) == CGAL::LEFT_TURN) {
            backward_shoot=true;
        }
        else{
            throw std::runtime_error("instance error: straight line consists of more than one edge");
        }

        bool forward_shoot;

        for(int i=0; i<n; i++){

            Point pt1(points[poly[i]].x, points[poly[i]].y);
            Point pt2(points[poly[((i+1) % n)]].x, points[poly[((i+1) % n)]].y);
            Point pt3(points[poly[((i+2) % n)]].x, points[poly[((i+2) % n)]].y);

            if (CGAL::orientation(pt1, pt2, pt3) == CGAL::RIGHT_TURN) {
                forward_shoot=false;

            } else if (CGAL::orientation(pt1, pt2, pt3) == CGAL::LEFT_TURN) {
                forward_shoot=true;
            } else{std::cerr<<"ERROR";}

            segments.emplace_back(Segment_w_info(Segment(pt1,pt2), true, seg_id++, poly_id, backward_shoot, forward_shoot));
            backward_shoot=forward_shoot;
        }
        poly_id++;
    }
}

void write_to_shp(std::vector<PWH> polys, std::string path) {

    std::string file = path + ".shp";
    //create handle
    SHPHandle shapefile = SHPCreate(file.c_str(), SHPT_POLYGON);
    DBFHandle dbfile = DBFCreate(file.c_str());

    //create field for poly ID
    int field_face_id = DBFAddField(dbfile, "ID", FTInteger, 5, 0);

    int f = 0;

    for (const auto& poly : polys) {
        // Count total number of vertices and parts
        int total_vertex_count = 0;
        int part_count = 1 + poly.number_of_holes(); // 1 for outer boundary + holes

        // Get the outer boundary and add its vertex count
        const Polygon_2& outer_boundary = poly.outer_boundary();
        total_vertex_count += outer_boundary.size();

        // Get hole vertex counts
        for (auto hole_it = poly.holes_begin(); hole_it != poly.holes_end(); ++hole_it) {
            total_vertex_count += hole_it->size();
        }


        //collect vertices
        double* x = new double[total_vertex_count];
        double* y = new double[total_vertex_count];

        // Array to store part start indices
        int* part_start_indices = new int[part_count];
        int vertex_index = 0;

        // Write outer boundary vertices
        part_start_indices[0] = 0; // Outer boundary starts at index 0
        for (auto vit = outer_boundary.vertices_begin(); vit != outer_boundary.vertices_end(); ++vit) {
            x[vertex_index] = to_double(vit->x());
            y[vertex_index] = to_double(vit->y());
            vertex_index++;
        }

        // Write each hole
        int part_index = 1;
        for (auto hole_it = poly.holes_begin(); hole_it != poly.holes_end(); ++hole_it) {
            part_start_indices[part_index++] = vertex_index;

            // Ensure the hole is CW
            Polygon_2 hole = *hole_it;
            if (hole.is_counterclockwise_oriented()) {
                hole.reverse_orientation();
            }

            for (auto vit = hole.vertices_begin(); vit != hole.vertices_end(); ++vit) {
                x[vertex_index] = to_double(vit->x());
                y[vertex_index] = to_double(vit->y());
                vertex_index++;
            }
        }

        //create polygon object
        SHPObject* shape = SHPCreateObject(SHPT_POLYGON, -1, part_count, part_start_indices, nullptr, total_vertex_count, x, y, nullptr, nullptr);

        //write shape into file
        int shape_id = SHPWriteObject(shapefile, -1, shape);

        //set field to face ID
        DBFWriteIntegerAttribute(dbfile, shape_id, field_face_id, f++);

        //memory management
        SHPDestroyObject(shape);
        delete[] x;
        delete[] y;
        delete[] part_start_indices;

    }

    DBFClose(dbfile);
    SHPClose(shapefile);
}


}

namespace GPKG {

inline std::string global_crs = "EPSG:25832";

void read(std::string filename, std::vector<Segment_w_info>& input_segments, Logger &logger) {
    auto polys = read_gpkg_to_pwh(filename);
    logger.add("Number of polygons in input", polys.size());
    pwh_to_swi(polys, input_segments);
}

std::vector<PWH> read_gpkg_to_pwh(std::string filename) {
    std::vector<PWH> polys;
    // Register all drivers
    GDALAllRegister();

    // Open the dataset
    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(
        filename.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));

    if (poDS == nullptr) {
        std::cerr << "Failed to open GPKG: " << filename << std::endl;
        return polys;
    }

    // Loop over layers
    for (int i = 0; i < poDS->GetLayerCount(); ++i) {
        OGRLayer* poLayer = poDS->GetLayer(i);
        if (poLayer == nullptr) continue;

        poLayer->ResetReading();
        OGRFeature* poFeature;
        while ((poFeature = poLayer->GetNextFeature()) != nullptr) {
            OGRGeometry* poGeometry = poFeature->GetGeometryRef();
            if (poGeometry != nullptr &&
                (wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon ||
                 wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon)) {

                std::vector<PWH> feature_polys;

                // Convert geometry to CGAL polygons
                if (wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon) {
                    PWH poly = OGRPolygonToCGAL(
                        dynamic_cast<OGRPolygon*>(poGeometry));
                    feature_polys.push_back(poly);
                } else if (wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon) {
                    OGRMultiPolygon* poMulti = dynamic_cast<OGRMultiPolygon*>(poGeometry);
                    for (auto&& poPart : *poMulti) {
                        PWH poly = OGRPolygonToCGAL(
                            dynamic_cast<OGRPolygon*>(poPart));
                        feature_polys.push_back(poly);
                    }
                }

                // Add polygons to output vector
                polys.insert(polys.end(), feature_polys.begin(), feature_polys.end());
                 }
            OGRFeature::DestroyFeature(poFeature);
        }

        // === get CRS and store it globally ===
        OGRSpatialReference* srs = poLayer->GetSpatialRef();
        if (srs != nullptr) {
            char* auth_name = nullptr;
            char* auth_code = nullptr;
            if (srs->AutoIdentifyEPSG() == OGRERR_NONE) {
                auth_name = const_cast<char*>(srs->GetAuthorityName(nullptr));
                auth_code = const_cast<char*>(srs->GetAuthorityCode(nullptr));
                if (auth_name && auth_code) {
                    global_crs = std::string(auth_name) + ":" + std::string(auth_code);
                }
            }
            if (global_crs.empty()) {
                char* wkt = nullptr;
                srs->exportToWkt(&wkt);
                if (wkt) {
                    global_crs = wkt;
                    CPLFree(wkt);
                }
            }
        }
    }

    GDALClose(poDS);
    return polys;
}

//writes polygons to geopackage
void write_to_gpkg(const std::vector<PWH>& polys, const std::string& path) {
    GDALAllRegister();

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GPKG");
    if (!driver) throw std::runtime_error("GPKG driver not available");

    std::string file = path + ".gpkg";
    GDALDataset* dataset = driver->Create(file.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!dataset) throw std::runtime_error("Failed to create GPKG file");

    OGRSpatialReference srs;
    if (srs.SetFromUserInput(global_crs.c_str()) != OGRERR_NONE) {
        GDALClose(dataset);
        throw std::runtime_error("Invalid CRS: " + global_crs);
    }
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    srs.AutoIdentifyEPSG();

    const char* lco[] = { "FID=id", "GEOMETRY_NAME=geom", nullptr };
    OGRLayer* layer = dataset->CreateLayer("", &srs, wkbPolygon, const_cast<char**>(lco));
    if (!layer) { GDALClose(dataset); throw std::runtime_error("Failed to create layer"); }

    for (const auto& poly : polys) {
        OGRPolygon ogr_poly;

        // Außenring
        OGRLinearRing outer;
        for (auto vit = poly.outer_boundary().vertices_begin(); vit != poly.outer_boundary().vertices_end(); ++vit) {
            outer.addPoint(to_double(vit->x()), to_double(vit->y()));
        }
        outer.closeRings();
        ogr_poly.addRing(&outer);

        // Löcher (Orientierung ggf. relativ zum Außenring anpassen)
        for (auto hole : poly.holes()) {
            OGRLinearRing inner;
            std::vector<Point> vs;
            for (auto v : hole.vertices()) vs.push_back(v);
            // Beispiel: Reverse nur, wenn nötig; hier nehmen wir umgekehrt wie Außenring an
            for (auto it = vs.rbegin(); it != vs.rend(); ++it)
                inner.addPoint(to_double(it->x()), to_double(it->y()));
            inner.closeRings();
            ogr_poly.addRing(&inner);
        }

        // Validieren/Reparieren
        OGRGeometry* to_write = nullptr;
        if (!ogr_poly.IsValid()) {
            OGRGeometry* fixed = ogr_poly.MakeValid();
            if (!fixed || !fixed->IsValid()) {
                if (fixed) OGRGeometryFactory::destroyGeometry(fixed);
                throw std::runtime_error("Invalid polygon geometry after MakeValid()");
            }
            to_write = fixed; // Besitz geht an Feature
        } else {
            to_write = ogr_poly.clone();
        }

        OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        f->SetGeometryDirectly(to_write);
        if (layer->CreateFeature(f) != OGRERR_NONE) {
            OGRFeature::DestroyFeature(f);
            GDALClose(dataset);
            throw std::runtime_error("Failed to create feature");
        }
        OGRFeature::DestroyFeature(f);
    }
    GDALClose(dataset);
}



//helper function to convert OGR polygons to CGAL polygons
PWH OGRPolygonToCGAL(OGRPolygon* poPolygon) {
    Polygon_2 outer;
    // Exterior ring
    OGRLinearRing* exterior = poPolygon->getExteriorRing();
    if (exterior) {
        int n = exterior->getNumPoints();
        if (n > 0) {
            // Add first point explicitly
            outer.push_back(Kernel::Point_2(exterior->getX(0), exterior->getY(0)));

            // Add subsequent points only if different from previous one
            for (int i = 1; i < n - 1; ++i) { // skip last point since it's a duplicate of the first
                double x = exterior->getX(i);
                double y = exterior->getY(i);
                Kernel::Point_2 curr(x, y);
                if (outer.size() == 0 || curr != outer[outer.size() - 1]) {
                    outer.push_back(curr);
                }
            }
        }
    }
    if (!outer.is_counterclockwise_oriented()) outer.reverse_orientation();

    std::vector<Polygon_2> holes;

    // Interior rings
    for (int i = 0; i < poPolygon->getNumInteriorRings(); ++i) {
        OGRLinearRing* ring = poPolygon->getInteriorRing(i);
        Polygon_2 hole;
        int n = ring->getNumPoints();
        if (n > 0) {
            hole.push_back(Kernel::Point_2(ring->getX(0), ring->getY(0)));

            for (int j = 1; j < n - 1; ++j) { // skip last point (duplicate of first)
                double x = ring->getX(j);
                double y = ring->getY(j);
                Kernel::Point_2 curr(x, y);
                if (hole.size() == 0 || curr != hole[hole.size() - 1]) {
                    hole.push_back(curr);
                }
            }
        }
        if (!hole.is_empty() && hole.is_simple()) {
            if (!hole.is_clockwise_oriented()) hole.reverse_orientation();
            holes.push_back(hole);
        }
    }

    //Error message on invalid polygons
    PWH poly(outer, holes.begin(), holes.end());

    return poly;
}


}


void pwh_to_swi(std::vector<PWH> &polygons, std::vector<Segment_w_info> &segments) {
    int seg_id=0;
    int poly_id=0;

    for (auto &poly:polygons) {
        Polygon_2 outer_boundary = poly.outer_boundary();

        int n = outer_boundary.size();
        if(n<2) continue;

        if(n==2){
            segments.emplace_back(Segment_w_info(Segment( outer_boundary[0],outer_boundary[1]), true, seg_id++, poly_id++, true, true));
            continue;
        }

        bool backward_shoot = (CGAL::orientation(outer_boundary[n - 1],outer_boundary[0],
            outer_boundary[1]) == CGAL::LEFT_TURN);

        bool forward_shoot;

        for(int i=0; i<n; i++) {

            forward_shoot = (CGAL::orientation(outer_boundary[i], outer_boundary[(i+1) % n], outer_boundary[(i+2) % n])
                == CGAL::LEFT_TURN);

            segments.emplace_back(Segment_w_info(Segment(outer_boundary[i],outer_boundary[(i+1) % n]), true, seg_id++, poly_id, backward_shoot, forward_shoot));

            backward_shoot=forward_shoot;

        }
        // Holes
        for (PWH::Hole_iterator hit = poly.holes_begin(); hit != poly.holes_end(); ++hit) {
            if (!hit->is_clockwise_oriented()) {
                hit->reverse_orientation();
                std::cout << "Had to reverse orientation of hole in polygon with ID: " << poly_id << std::endl;
            }

            const auto &hole = *hit;
            if (hole.size() < 3) continue;

            int m = hole.size();

            backward_shoot = (CGAL::orientation(hole[m - 1],hole[0],
                hole[1]) == CGAL::LEFT_TURN);

            for(int j=0; j<m; j++) {

                forward_shoot = (CGAL::orientation(hole[j], hole[(j+1) % m], hole[(j+2) % m])
                    == CGAL::LEFT_TURN);

                segments.emplace_back(Segment_w_info(Segment(hole[j],hole[(j+1) % m]), false, seg_id++,
                    poly_id, backward_shoot, forward_shoot));
                backward_shoot=forward_shoot;
            }
        }
        poly_id++;
    }
}
namespace SVG {
    void segments_to_svg(const std::vector<Segment>& segments, const std::string& filename) {

        double minX = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        double maxY = std::numeric_limits<double>::lowest();

        for (const auto& s : segments) {
            minX = std::min({minX, CGAL::to_double(s.source().x()), CGAL::to_double(s.target().x())});
            maxX = std::max({maxX, CGAL::to_double(s.source().x()), CGAL::to_double(s.target().x())});
            minY = std::min({minY, CGAL::to_double(s.source().y()), CGAL::to_double(s.target().y())});
            maxY = std::max({maxY, CGAL::to_double(s.source().y()), CGAL::to_double(s.target().y())});
        }

        double padding = 10.0; // Adjust based on your coordinate range
        minX -= padding;
        maxX += padding;
        minY -= padding;
        maxY += padding;

        double width = maxX - minX;
        double height = maxY - minY;

        std::ofstream svg_file(filename);
        svg_file << std::fixed << std::setprecision(10); // HIGH precision output
        double pixel_width = 1000.0;
        double aspect_ratio = height / width;
        double pixel_height = pixel_width * aspect_ratio;
        double scale = pixel_width / width; // pixel per unit
        double stroke_width_units = 1.0 / scale; // 1 screen pixel = X world units

        svg_file << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
                 << "viewBox=\"" << minX << " " << minY << " " << width << " " << height << "\" "
                 << "width=\"" << pixel_width << "\" height=\"" << pixel_height << "\" "
                 << "preserveAspectRatio=\"xMidYMid meet\" "
                 << "fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width_units << "\">\n";


        for (const auto& s : segments) {
            svg_file << "<line x1=\"" << s.source().x() << "\" y1=\"" << (maxY -(s.source().y()- minY))
                     << "\" x2=\"" << s.target().x() << "\" y2=\"" << (maxY - (s.target().y() - minY))
                     << "\" />\n";
        }

        svg_file << "</svg>\n";

        svg_file.close();
    }

    void polygons_to_svg(const std::vector<PWH>& polygons, const std::string& filename) {

        double minX = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        double maxY = std::numeric_limits<double>::lowest();

        // Compute bounding box over all outer and hole polygons
        for (const auto& pwh : polygons) {
            for (const auto& p : pwh.outer_boundary()) {
                double x = CGAL::to_double(p.x());
                double y = CGAL::to_double(p.y());
                minX = std::min(minX, x);
                maxX = std::max(maxX, x);
                minY = std::min(minY, y);
                maxY = std::max(maxY, y);
            }
            for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
                for (const auto& p : *hit) {
                    double x = CGAL::to_double(p.x());
                    double y = CGAL::to_double(p.y());
                    minX = std::min(minX, x);
                    maxX = std::max(maxX, x);
                    minY = std::min(minY, y);
                    maxY = std::max(maxY, y);
                }
            }
        }

        double padding = 10.0;
        minX -= padding;
        maxX += padding;
        minY -= padding;
        maxY += padding;

        double width = maxX - minX;
        double height = maxY - minY;

        std::ofstream svg_file(filename);
        svg_file << std::fixed << std::setprecision(10);

        double pixel_width = 1000.0;
        double aspect_ratio = height / width;
        double pixel_height = pixel_width * aspect_ratio;
        double scale = pixel_width / width;
        double stroke_width_units = 1.0 / scale;

        svg_file << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
                 << "viewBox=\"" << minX << " " << minY << " " << width << " " << height << "\" "
                 << "width=\"" << pixel_width << "\" height=\"" << pixel_height << "\" "
                 << "preserveAspectRatio=\"xMidYMid meet\" "
                 << "fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width_units << "\">\n";

        for (const auto& pwh : polygons) {
            // Outer boundary
            svg_file << "<polygon points=\"";
            for (const auto& p : pwh.outer_boundary()) {
                double x = CGAL::to_double(p.x());
                double y = CGAL::to_double(p.y());
                svg_file << x << "," << (maxY - (y - minY)) << " ";
            }
            svg_file << "\" fill=\"lightgray\" stroke=\"black\"/>\n";

            // Holes
            for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
                svg_file << "<polygon points=\"";
                for (const auto& p : *hit) {
                    double x = CGAL::to_double(p.x());
                    double y = CGAL::to_double(p.y());
                    svg_file << x << "," << (maxY - (y - minY)) << " ";
                }
                svg_file << "\" fill=\"white\" stroke=\"black\"/>\n";
            }
        }

        svg_file << "</svg>\n";
        svg_file.close();
    }
}// namespace SVG

void all_polygons_in_solution(std::vector<PWH> &    all_polys_in_sol, std::vector<bool> max_flow_solution, const Arrangement& arr) {
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (!max_flow_solution[fit->data().id]) {
            continue;
        }
        Polygon_2 poly;
        // outer boundary exists
        if (!fit->has_outer_ccb()) {
            continue;
        }
        auto circ = fit->outer_ccb();
        auto curr = circ;
        do {
            const auto& source = curr->source()->point();
            poly.push_back(source);
            ++curr;
        } while (curr != circ);
        PWH pwh(poly);
        all_polys_in_sol.push_back(pwh);
    }
}

void arrangement_as_polys(std::vector<PWH> &all_polys, const Arrangement& arr) {
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        Polygon_2 poly;
        // outer boundary exists
        if (!fit->has_outer_ccb()) {
            continue;
        }
        auto circ = fit->outer_ccb();
        auto curr = circ;
        do {
            const auto& source = curr->source()->point();
            poly.push_back(source);
            ++curr;
        } while (curr != circ);
        PWH pwh(poly);
        all_polys.push_back(pwh);
    }
}



const std::pair <std::vector<Polygon_2>, std::vector<Polygon_2>>  combine_polygons(const std::vector<bool>
&in_solution, Arrangement &arr) {
    std::cout<<"GROUPSSIZE"<<in_solution.size()<<std::endl;
    std::vector<Polygon_2> outer_boundaries;
    std::vector<Polygon_2> holes;
    for (Arrangement::Halfedge_iterator eit = arr.halfedges_begin(); eit != arr.halfedges_end(); ++eit) {
        auto eitc = eit;
        if (!eitc->data().visited && in_solution[eitc->face()->data().id] && (!in_solution[eitc->twin()->face()->data().id]||eitc->twin()->face()->is_unbounded())) {
            Polygon_2 poly = get_contiguous_boundary(eitc, in_solution);
            assert(poly.is_simple());
            if (poly.is_counterclockwise_oriented()) {
                outer_boundaries.emplace_back(poly);
            } else {
                holes.emplace_back(poly);
            }
        }

    }
    auto holes_and_outer = std::make_pair(outer_boundaries, holes);
    return holes_and_outer;

}


const Polygon_2 get_contiguous_boundary(Arrangement::Halfedge_handle &edge, const std::vector<bool> &in_solution) {

    assert(in_solution[edge->face()->data().id] && (!in_solution[edge->twin()->face()->data().id])||edge->twin()->face()->is_unbounded());
    Polygon_2 polygon;
    do {
        edge->data().visited = true;
        if (!in_solution[edge->twin()->face()->data().id] || edge->twin()->face()->is_unbounded()) {
            polygon.push_back(edge->source()->point());
            if (edge->target()->point() == *(polygon.vertices_begin())) {
                return polygon; // finished
            }
            edge = edge->next();
        }
        else {
            edge = edge->twin()->next();
        }
    }while (true);

}

const std::map<int, std::vector<Polygon_2>> locate_holes(const std::vector<Polygon_2> &outer_boundaries,
    const std::vector<Polygon_2> &holes) {

    std::vector<CGAL::Bbox_2> outer_bboxes;
    for (const auto& poly : outer_boundaries) {
        outer_bboxes.emplace_back(poly.bbox());
    }

    std::map<int, std::vector<Polygon_2>> outer_to_holes;
    for (int i = 0; i < holes.size(); ++i) {
        const Polygon_2& hole = holes[i];
        Point test_point = *hole.vertices_begin();
        CGAL::Bbox_2 point_bbox = test_point.bbox();

        for (int j = 0; j < outer_boundaries.size(); ++j) {
            if (!CGAL::do_overlap(point_bbox, outer_bboxes[j]))
                continue;

            CGAL::Bounded_side side = CGAL::bounded_side_2(
                outer_boundaries[j].vertices_begin(),
                outer_boundaries[j].vertices_end(),
                test_point,
                Kernel());

            if (side == CGAL::ON_BOUNDED_SIDE) {
                outer_to_holes[j].emplace_back(holes[i]);  // hole i is inside outer j
                break;
            }
        }
    }
    return outer_to_holes;
}

const std::vector<PWH> create_polygons_with_holes(const std::vector<Polygon_2> &outer_boundaries, const std::vector<Polygon_2> &holes) {
    auto outer_to_holes = locate_holes(outer_boundaries, holes);

    std::vector<PWH> polygons;

    std::vector<bool> already_added (outer_boundaries.size(), false);
    for (const auto &outer_pair : outer_to_holes) {
        int compare = polygons.size();
        PWH pwh(outer_boundaries[outer_pair.first], outer_pair.second.begin(), outer_pair.second.end());
        polygons.emplace_back(pwh);
        already_added[outer_pair.first] = true;
    }
    for (int outer_boundary = 0; outer_boundary < outer_boundaries.size(); ++outer_boundary) {
        if (!already_added[outer_boundary]) {
            int compare = polygons.size();

            PWH pwh(outer_boundaries[outer_boundary]);
            polygons.emplace_back(pwh);

        }
    }
    return polygons;
}

const std::vector<PWH> cgal_combines(const std::vector<Polygon_2> &polygons){

    CGAL::Polygon_set_2<Kernel> pset;
    pset.join(polygons.begin(), polygons.end());
    std::vector<PWH> out;
    pset.polygons_with_holes(std::back_inserter(out));
    return out;
}

}



