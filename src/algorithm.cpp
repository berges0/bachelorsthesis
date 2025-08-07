#include "algorithm.hpp"


void algorithm::run(std::string input_path, double alpha) {

    // Load instance using the library
    auto data = SHPLoader::ReadShapeFileToPoint2D(input_path);
    std::cout << "Successfully loaded instance from: " << input_path << std::endl;
    std::cout << "Points: " << data.first.size() << std::endl;
    std::cout << "Polygons: " << data.second.size() << std::endl;
    std::pair<std::vector<Segment_w_info>, std::pair<Point, Point>> transformed_to_segments;

    bool debugging = false;
    if (debugging) {
        std::vector<SHPLoader::Point2D> deb_points = {
            SHPLoader::Point2D(0, 0), SHPLoader::Point2D(4, 0), SHPLoader::Point2D(4, 4), SHPLoader::Point2D(2.5,4), SHPLoader::Point2D(3,1)
            , SHPLoader::Point2D(1,1), SHPLoader::Point2D(1.5,4), SHPLoader::Point2D(0,4),


            SHPLoader::Point2D(-3, 5), SHPLoader::Point2D(10, 5), SHPLoader::Point2D(10,7), SHPLoader::Point2D(-3,7),

            SHPLoader::Point2D(5, 2), SHPLoader::Point2D(9, 2), SHPLoader::Point2D(9, 8), SHPLoader::Point2D(5, 8),  SHPLoader::Point2D(5, 7),  SHPLoader::Point2D(6, 7),  SHPLoader::Point2D(6, 5),  SHPLoader::Point2D(5, 5),

            // Rechteck – Gebäude 1
            SHPLoader::Point2D(0, 0), SHPLoader::Point2D(0, 4), SHPLoader::Point2D(4, 4), SHPLoader::Point2D(4, 0),

            // Fünfeck – Gebäude 2
            SHPLoader::Point2D(10, 0), SHPLoader::Point2D(12, 4), SHPLoader::Point2D(14, 4), SHPLoader::Point2D(16, 0), SHPLoader::Point2D(13, -2),


            // Sechseck – Gebäude 3
            SHPLoader::Point2D(20, 0), SHPLoader::Point2D(22, 2), SHPLoader::Point2D(24, 2), SHPLoader::Point2D(26, 0), SHPLoader::Point2D(24, -2), SHPLoader::Point2D(22, -2),

            // L-förmiges Polygon – Gebäude 4
            SHPLoader::Point2D(0, 10), SHPLoader::Point2D(0, 14), SHPLoader::Point2D(2, 14), SHPLoader::Point2D(2, 12), SHPLoader::Point2D(4, 12), SHPLoader::Point2D(4, 10),

            // Trapez – Gebäude 5
            SHPLoader::Point2D(10, 10), SHPLoader::Point2D(12, 14), SHPLoader::Point2D(16, 14), SHPLoader::Point2D(18, 10),


            // Rechteck – Gebäude 6
            SHPLoader::Point2D(20, 10), SHPLoader::Point2D(20, 14), SHPLoader::Point2D(25, 14), SHPLoader::Point2D(25, 10),

            // Unregelmäßiges Fünfeck – Gebäude 7
            SHPLoader::Point2D(0, 20), SHPLoader::Point2D(1, 24), SHPLoader::Point2D(3, 25), SHPLoader::Point2D(5, 23), SHPLoader::Point2D(4, 20),

            // Sechseck – Gebäude 8
            SHPLoader::Point2D(10, 20), SHPLoader::Point2D(11, 22), SHPLoader::Point2D(13, 22), SHPLoader::Point2D(14, 20), SHPLoader::Point2D(13, 18), SHPLoader::Point2D(11, 18),

            // Rechteck – Gebäude 9
            SHPLoader::Point2D(20, 20), SHPLoader::Point2D(20, 24), SHPLoader::Point2D(24, 24), SHPLoader::Point2D(24, 20),

            // Dreieck – Gebäude 10
            SHPLoader::Point2D(30, 0), SHPLoader::Point2D(32, 4), SHPLoader::Point2D(34, 0),

            // Fünfeck – Gebäude 11
            SHPLoader::Point2D(30, 10), SHPLoader::Point2D(32, 13), SHPLoader::Point2D(35, 13), SHPLoader::Point2D(37, 10), SHPLoader::Point2D(34, 8),

            // Rechteck – Gebäude 12
            SHPLoader::Point2D(30, 20), SHPLoader::Point2D(30, 25), SHPLoader::Point2D(34, 25), SHPLoader::Point2D(34, 20),

            // Sechseck – Gebäude 13
            SHPLoader::Point2D(0, 30), SHPLoader::Point2D(2, 32), SHPLoader::Point2D(4, 32), SHPLoader::Point2D(6, 30), SHPLoader::Point2D(4, 28), SHPLoader::Point2D(2, 28),

            // Trapez – Gebäude 14
            SHPLoader::Point2D(10, 30), SHPLoader::Point2D(12, 34), SHPLoader::Point2D(16, 34), SHPLoader::Point2D(18, 30),

            // Rechteck – Gebäude 15
            SHPLoader::Point2D(20, 30), SHPLoader::Point2D(20, 34), SHPLoader::Point2D(24, 34), SHPLoader::Point2D(24, 30),

            // Fünfeck – Gebäude 16
            SHPLoader::Point2D(30, 30), SHPLoader::Point2D(32, 34), SHPLoader::Point2D(35, 34), SHPLoader::Point2D(37, 30), SHPLoader::Point2D(33.5, 28),

            // Rechteck – Gebäude 17
            SHPLoader::Point2D(40, 0), SHPLoader::Point2D(40, 5), SHPLoader::Point2D(45, 5), SHPLoader::Point2D(45, 0),

            // Sechseck – Gebäude 18
            SHPLoader::Point2D(40, 10), SHPLoader::Point2D(42, 12), SHPLoader::Point2D(44, 12), SHPLoader::Point2D(46, 10), SHPLoader::Point2D(44, 8), SHPLoader::Point2D(42, 8),

            // Trapez – Gebäude 19
            SHPLoader::Point2D(40, 20), SHPLoader::Point2D(42, 23), SHPLoader::Point2D(46, 23), SHPLoader::Point2D(48, 20),

            // Dreieck – Gebäude 20
            SHPLoader::Point2D(40, 30), SHPLoader::Point2D(42, 34), SHPLoader::Point2D(44, 30)
        };

        std::vector<std::vector<int>> deb_polygons = {
            {0,1,2,3,4,5,6,7},
            {8,9,10,11}
            /*{4,5,6,7,8,9,11}
            {0, 3, 2, 1},                // Gebäude 1 – Rechteck
            {4, 8, 7, 6, 5},             // Gebäude 2 – Fünfeck
            {9, 14, 13, 12, 11, 10},     // Gebäude 3 – Sechseck
            {15, 20, 19, 18, 17, 16},    // Gebäude 4 – L-Form
            {21, 24, 23, 22},            // Gebäude 5 – Trapez
            {25, 28, 27, 26},            // Gebäude 6 – Rechteck
            {29, 33, 32, 31, 30},        // Gebäude 7 – Unreg. Fünfeck
            {34, 39, 38, 37, 36, 35},    // Gebäude 8 – Sechseck
            {40, 43, 42, 41},            // Gebäude 9 – Rechteck
            {44, 46, 45},                // Gebäude 10 – Dreieck
            {47, 51, 50, 49, 48},        // Gebäude 11 – Fünfeck
            {52, 55, 54, 53},            // Gebäude 12 – Rechteck
            {56, 61, 60, 59, 58, 57},    // Gebäude 13 – Sechseck
            {62, 65, 64, 63},            // Gebäude 14 – Trapez
            {66, 69, 68, 67},            // Gebäude 15 – Rechteck
            {70, 74, 73, 72, 71},        // Gebäude 16 – Fünfeck
            {75, 78, 77, 76},            // Gebäude 17 – Rechteck
            {79, 84, 83, 82, 81, 80},    // Gebäude 18 – Sechseck
            {85, 88, 87, 86},            // Gebäude 19 – Trapez
            {89, 91, 90}                 // Gebäude 20 – Dreieck*/
        };
        transformed_to_segments = build_segment_vector(deb_points, deb_polygons);
    }
    else {
        transformed_to_segments = build_segment_vector(data.first, data.second);
    }
    //auto transformed_to_segments = build_segment_vector(points, polygons);
    std::vector<Segment_w_info> segments_w_info = transformed_to_segments.first;
    add_box(segments_w_info, transformed_to_segments.second);
    std::vector<Segment> segments= segs_wo_info(segments_w_info);
    segments_to_svg(segments,"original.svg");
    auto extended_segs = ray_shoot_intersection(segments_w_info);
    std::cout<<"Extended Segments size: "<<extended_segs.size()<<std::endl;
    segments_to_svg(segs_wo_info(extended_segs),"extended.svg");
    
    Arrangement arr = build_arrangement(extended_segs);
    auto graph = build_graph(arr, alpha);
    auto max_sol = max_flow(graph,arr);
    auto output_data = combined_output_polygons(max_sol, arr);
/*
    Polygon_2 outer;
    outer.push_back(Point(0, 0));
    outer.push_back(Point(10, 0));
    outer.push_back(Point(10, 10));
    outer.push_back(Point(0, 10));

    // Define hole (smaller square from (3,3) to (7,7))
    Polygon_2 hole;
    hole.push_back(Point(3, 3));
    hole.push_back(Point(7, 3));
    hole.push_back(Point(7, 7));
    hole.push_back(Point(3, 7));

    // Ensure correct orientation
    if (!outer.is_counterclockwise_oriented())
        outer.reverse_orientation();
    if (!hole.is_clockwise_oriented())
        hole.reverse_orientation();

    // Construct polygon with hole
    Polygon_with_holes_2 poly_with_hole(outer);
    poly_with_hole.add_hole(hole);
    std::vector<Polygon_with_holes_2> poly_with_holes = {poly_with_hole};
    polygons_to_svg(poly_with_holes,"test_poly_with_holes.svg");
*/
    polygonns_to_svg(output_polygons(max_sol, arr), "output_non_combined.svg");
    //auto output_data = output_polygons(max_sol, arr);
    polygons_to_svg(output_data, "output_polys.svg");
    write_to_shp(output_data, "/home/samuel-berge/Work/bachelorsthesis/build/shp/shp_output.shp");
    //write_output_polygons(output_data, out_path);
    std::cout << "Arrangement has "
          << arr.number_of_vertices() << " vertices, "
          << arr.number_of_edges() << " edges, "
          << arr.number_of_faces() << " faces." << std::endl;
}