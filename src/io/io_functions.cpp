//
// Created by samuel-berge on 8/11/25.
//
#include "io/io_functions.hpp"
//
// Created by samuel-berge on 8/10/25.
//

namespace IO_FUNCTIONS {

//CONTROVERSION FROM SHAPEFILEAGGREGATIONLOADER TO CGAL DATASTRUCTURE

std::vector<Segment_w_info> read_in_segment(const std::vector<SHPLoader::Point2D>& points,
    const std::vector<std::vector<int>>& polygons) {

    std::vector<Segment_w_info> segments;

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
    return segments;
}


//SVG FILES

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

void polygons_to_svg(const std::vector<Polygon_with_holes_2>& polygons, const std::string& filename) {

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

//SHP FILES

void write_to_shp(const std::vector<Polygon_with_holes_2>& polygons, const std::string& filename) {

    std::vector<std::pair<double, double>> points;
    std::vector<std::vector<int>> polys;

    std::vector<int> polygon_indices;

    for (const auto& pwh : polygons) {
        // Outer boundary
        polygon_indices.clear();
        for (const auto& p : pwh.outer_boundary()) {
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            points.emplace_back(x, y);
            polygon_indices.emplace_back(points.size() - 1);
        }
        polys.push_back(polygon_indices);

        // Holes
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
            polygon_indices.clear();
            for (const auto& p : *hit) {
                double x = CGAL::to_double(p.x());
                double y = CGAL::to_double(p.y());
                points.emplace_back(x, y);
                polygon_indices.emplace_back(points.size() - 1);
            }
            polys.push_back(polygon_indices);
        }
    }

    std::cout << "THERE ARE " << points.size() << " points" << std::endl;
    std::cout << "AND " << polys.size() << " polygons (including holes)" << std::endl;

    SHPLoader::writeToShapeFile(std::make_pair(points, polys), filename);
}

//CONVERT SOLUTION POLYGONS TO CONTIGUOUS POLYGONS


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

const std::vector<Polygon_with_holes_2> create_polygons_with_holes(const std::vector<Polygon_2> &outer_boundaries, const std::vector<Polygon_2> &holes) {
    auto outer_to_holes = locate_holes(outer_boundaries, holes);

    std::vector<Polygon_with_holes_2> polygons;

    std::vector<bool> already_added (outer_boundaries.size(), false);
    for (const auto &outer_pair : outer_to_holes) {
        int compare = polygons.size();
        Polygon_with_holes_2 pwh(outer_boundaries[outer_pair.first], outer_pair.second.begin(), outer_pair.second.end());
        polygons.emplace_back(pwh);
        already_added[outer_pair.first] = true;
    }
    for (int outer_boundary = 0; outer_boundary < outer_boundaries.size(); ++outer_boundary) {
        if (!already_added[outer_boundary]) {
            int compare = polygons.size();

            Polygon_with_holes_2 pwh(outer_boundaries[outer_boundary]);
            polygons.emplace_back(pwh);

        }
    }
    return polygons;
}

const std::vector<Polygon_with_holes_2> cgal_combines(const std::vector<Polygon_2> &polygons){

    CGAL::Polygon_set_2<Kernel> pset;
    pset.join(polygons.begin(), polygons.end());
    std::vector<Polygon_with_holes_2> out;
    pset.polygons_with_holes(std::back_inserter(out));
    return out;
}
void writeToShapeFile(std::vector<Polygon_with_holes_2> polys, std::string path) {

    //create handle
    SHPHandle shapefile = SHPCreate(path.c_str(), SHPT_POLYGON);
    DBFHandle dbfile = DBFCreate(path.c_str());

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


