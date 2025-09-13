//
// Created by samuel-berge on 8/20/25.
//
#include "core/algorithm.hpp"


namespace ALGORITHM{

void read_in (std::vector<Segment_w_info> &input_segments, const std::string &input_filename, Logger &logger) {

    auto pos = input_filename.find_last_of('.');
    if (pos == std::string::npos) {
        throw std::runtime_error("input file doesnt have an extension");
    };

    std::string ext = input_filename.substr(pos + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext == "shp") {
        IO_FUNCTIONS::SHP::read(input_filename, input_segments, logger);
    }
    else if (ext == "gpkg") {
        IO_FUNCTIONS::GPKG::read(input_filename, input_segments, logger);
        //std::filesystem::copy_file(input_filename, logger.out_dir_stem()+"_input.gpkg");
    }
    logger.add("Number of segments in input", input_segments.size());
    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(input_segments), "input.svg" );
}

void aggregate(std::vector<PWH> &output_data, const Arrangement &arr, Logger &logger) {

    double alpha = logger.alpha();

    logger.start_operation();
    Graph graph = GRAPH::build_graph(arr, alpha, logger);
    logger.end_operation("Building Graph (milliseconds)");
    logger.add("Number of edges in graph", std::get<0>(graph).size());


    logger.start_operation();
    std::vector<bool> max_flow_solution = MAX_FLOW::max_flow(graph, arr);
    logger.end_operation("Running Max Flow (milliseconds)");
    int nr_faces_solution = 0;
    for (int i = 0; i < max_flow_solution.size(); i++) {if (max_flow_solution[i]) nr_faces_solution++;}
    logger.add("Number of faces solution", nr_faces_solution);

    logger.start_operation();
    //IO_FUNCTIONS::all_polygons_in_solution(output_data, max_flow_solution, arr);
    auto holes_and_outer = IO_FUNCTIONS::combine_polygons(max_flow_solution, arr);
    output_data = IO_FUNCTIONS::create_polygons_with_holes(holes_and_outer.first, holes_and_outer.second);
    logger.end_operation("Combining faces of solution (milliseconds)");
    logger.add("Number of polygons in solution", holes_and_outer.first.size());
    logger.add("Number of holes in solution", holes_and_outer.second.size());
}

void run_standard(const std::string &input_filename, const std::string &output_filename, double alpha, Logger &logger) {

    std::vector<Segment_w_info> input_segments;

    read_in(input_segments, input_filename, logger);

    EDGE_EXTENSION::add_outer_box(input_segments,0.01);

    logger.start_operation();
    std::vector<Segment_w_info> extended_segments = EDGE_EXTENSION::STANDARD::extension(input_segments);
    logger.end_operation("Extending edges (milliseconds) ");


    logger.start_operation();
    Arrangement arr = ARRANGEMENT::build_arrangement(extended_segments, logger);
    logger.end_operation("Building Arragnement (milliseconds)");
    logger.add("Number of edges in arrangement", arr.number_of_edges());
    logger.add("Number of halfedges in arrangement", arr.number_of_halfedges());
    logger.add("Number of faces in arrangement", arr.number_of_faces());

    logger.stop_time();
    std::vector<PWH> output_data;
    std::vector<double> alpha_values= {1.0, 0.5, 0.25, 0.15, 0.1 ,0.05, 0.035, 0.025, 0.02,  0.01 ,0.001 ,0.005 ,0.00025,
        0.0001 , 0.00001, 0.0};


    std::vector<Logger> loggers;
    for (double alpha : alpha_values) {

        Logger logger_sub(logger,alpha);
        logger_sub.add("Alpha", logger_sub.alpha());

        output_data.clear();
        logger_sub.start();

        aggregate(output_data, arr, logger_sub);

        //IO_FUNCTIONS::writeToShapeFile(output_data, output_filename);

        IO_FUNCTIONS::SHP::write_to_shp(output_data, logger_sub.out_dir_stem() + "_solution");

        //std::string command = "qgis " + logger_sub.out_dir_stem() + "_solution.gpkg " + input_filename + " &";
        //int status = std::system(command.c_str());

        logger_sub.end_sub_log();
    }
}

void run_limited(const std::string &input_filename, const std::string &output_filename, double alpha, double threshold, double th_scale,
    int th_variant, Logger &logger) {

    std::vector<Segment_w_info> input_segments;

    read_in(input_segments, input_filename, logger);

    EDGE_EXTENSION::add_outer_box(input_segments,0.01);

    logger.start_operation();
    std::vector<Segment_w_info> extended_segments = EDGE_EXTENSION::LIMITED::extension(input_segments, threshold,
        th_scale, th_variant);
    logger.end_operation("Extending edges (milliseconds) ");

    logger.start_operation();
    Arrangement arr = ARRANGEMENT::build_arrangement(extended_segments, logger);
    logger.end_operation("Building Arragnement (milliseconds)");
    logger.add("Number of edges in arrangement", arr.number_of_edges());
    logger.add("Number of halfedges in arrangement", arr.number_of_halfedges());
    logger.add("Number of faces in arrangement", arr.number_of_faces());

    logger.stop_time();
    std::vector<PWH> output_data;
    std::vector<double> alpha_values= {1.0, 0.5, 0.25, 0.15, 0.1 ,0.05, 0.035, 0.025, 0.02,  0.01 ,0.001 ,0.005 ,0.00025,
        0.0001 , 0.00001, 0.0};

    std::vector<Logger> loggers;

    for (double alpha : alpha_values) {
        Logger logger_sub(logger,alpha);

        logger_sub.add("Alpha", logger_sub.alpha());

        output_data.clear();
        logger_sub.start();


        aggregate(output_data, arr, logger_sub);

        //IO_FUNCTIONS::writeToShapeFile(output_data, output_filename);

        IO_FUNCTIONS::SHP::write_to_shp(output_data, logger_sub.out_dir_stem() + "_solution");

        //std::string command = "qgis " + logger_sub.out_dir_stem() + "_solution.gpkg " + input_filename + " &";
        //int status = std::system(command.c_str());

        logger_sub.end_sub_log();
    }
}

void run_edge_relink(const std::string &input_filename, const std::string &output_filename, double alpha, double threshold, double th_scale,
    int th_variant, double degree, double distance, std::string subversion, Logger &logger) {

    std::vector<Segment_w_info> input_segments;

    read_in(input_segments, input_filename, logger);

    auto polygonswh = IO_FUNCTIONS::GPKG::read_gpkg_to_pwh(input_filename);

    RTree rtree = SUBSTITUTE_EDGES::build_r_tree(polygonswh);

    EDGE_EXTENSION::add_outer_box(input_segments,0.01);

    logger.start_operation();
    std::vector<Segment_w_info> extended_segments;
    if (subversion == "0") {
        extended_segments = EDGE_EXTENSION::STANDARD::extension(input_segments);
    }
    else if (subversion == "1") {
        extended_segments = EDGE_EXTENSION::LIMITED::extension(input_segments, threshold, th_scale, th_variant );
    }
    else {
        throw std::runtime_error("subversion not recognized");
    }
    logger.end_operation("Extending edges (milliseconds) ");
    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(extended_segments), "after_extension.svg");


    std::cout<<"Number of segments after extension: "<<extended_segments.size()<<std::endl;


    auto spatially_close = PRE_PROCESS::group_by_degree_and_closeness(extended_segments, degree,distance);

    auto just_seg = EDGE_EXTENSION::filter_segments(input_segments);
    Tree tree(just_seg.begin(), just_seg.end());
    tree.build();
    tree.accelerate_distance_queries();

    SUBSTITUTE_EDGES::relink_edges(spatially_close, tree);

    auto relinked_segments = spatially_close[0];

    if (subversion=="1") {
        SUBSTITUTE_EDGES::post_prune(relinked_segments);
    }
    std::cout<<"Number of segments after relinking: "<<relinked_segments.size()<<std::endl;

    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(relinked_segments), "after_relink.svg");


    std::vector<PWH> output_data;


    logger.start_operation();
    Arrangement arr = ARRANGEMENT::build_arrangement_relinked(relinked_segments, rtree, polygonswh, logger);
    logger.end_operation("Building Arragnement (milliseconds)");
    logger.add("Number of edges in arrangement", arr.number_of_edges());
    logger.add("Number of halfedges in arrangement", arr.number_of_halfedges());
    logger.add("Number of faces in arrangement", arr.number_of_faces());

    logger.stop_time();
    std::vector<double> alpha_values= {1.0, 0.5, 0.25, 0.15, 0.1 ,0.05, 0.035, 0.025, 0.02,  0.01 ,0.001 ,0.005 ,0.00025,
        0.0001 , 0.00001, 0.0};

    std::vector<Logger> loggers;

    for (double alpha : alpha_values) {
        Logger logger_sub(logger,alpha);
        logger_sub.add("Alpha", logger_sub.alpha());

        output_data.clear();
        logger_sub.start();

        aggregate(output_data, arr, logger_sub);

        //IO_FUNCTIONS::writeToShapeFile(output_data, output_filename);

        IO_FUNCTIONS::SHP::write_to_shp(output_data, logger_sub.out_dir_stem() + "_solution");

        //std::string command = "qgis " + logger_sub.out_dir_stem() + "_solution.gpkg " + input_filename + " &";
        //int status = std::system(command.c_str());

        logger_sub.end_sub_log();
    }
}


void run_outer_endpoints(const std::string &input_filename, const std::string &output_filename, double alpha, double threshold,
    double th_scale, int th_variant, double degree, double distance, std::string subversion, Logger &logger) {
    std::vector<Segment_w_info> input_segments;

    read_in(input_segments, input_filename, logger);

    auto polygonswh = IO_FUNCTIONS::GPKG::read_gpkg_to_pwh(input_filename);

    RTree rtree = SUBSTITUTE_EDGES::build_r_tree(polygonswh);

    EDGE_EXTENSION::add_outer_box(input_segments,0.01);

    logger.start_operation();
    std::vector<Segment_w_info> extended_segments;
    if (subversion == "0") {
        extended_segments = EDGE_EXTENSION::STANDARD::extension(input_segments);
    }
    else if (subversion == "1") {
        extended_segments = EDGE_EXTENSION::LIMITED::extension(input_segments, threshold, th_scale, th_variant);
    }
    else {
        throw std::runtime_error("subversion not recognized");
    }
    logger.end_operation("Extending edges (milliseconds) ");
    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(extended_segments), "after_extension.svg");


    std::cout<<"Number of segments after extension: "<<extended_segments.size()<<std::endl;

    auto spatially_close = PRE_PROCESS::group_by_degree_and_closeness(extended_segments, degree,distance);

    SUBSTITUTE_EDGES::connect_outer_points(spatially_close);

    auto connected_outer_point = spatially_close[0];

    if (subversion=="1") {
        SUBSTITUTE_EDGES::post_prune(connected_outer_point);
    }

    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(connected_outer_point), "after_connect_outer.svg");

    std::cout<<"Number of segments after relinking: "<<connected_outer_point.size()<<std::endl;


    logger.start_operation();
    Arrangement arr = ARRANGEMENT::build_arrangement_relinked(connected_outer_point, rtree, polygonswh, logger);
    logger.end_operation("Building Arragnement (milliseconds)");
    logger.add("Number of edges in arrangement", arr.number_of_edges());
    logger.add("Number of halfedges in arrangement", arr.number_of_halfedges());
    logger.add("Number of faces in arrangement", arr.number_of_faces());

    logger.stop_time();
    std::vector<PWH> output_data;
    std::vector<double> alpha_values= {1.0, 0.5, 0.25, 0.15, 0.1 ,0.05, 0.035, 0.025, 0.02,  0.01 ,0.001 ,0.005 ,0.00025,
        0.0001 , 0.00001, 0.0};

    std::vector<Logger> loggers;

    for (double alpha : alpha_values) {
        Logger logger_sub(logger,alpha);

        logger_sub.add("Alpha", logger_sub.alpha());

        output_data.clear();
        logger_sub.start();

        aggregate(output_data, arr, logger_sub);

        //IO_FUNCTIONS::writeToShapeFile(output_data, output_filename);

        IO_FUNCTIONS::SHP::write_to_shp(output_data, logger_sub.out_dir_stem() + "_solution");

        std::string command = "qgis " + logger_sub.out_dir_stem() + "_solution.gpkg " + input_filename + " &";
        int status = std::system(command.c_str());

        logger_sub.end_sub_log();
    }
}

void run_subdivision(const std::string &input_filename, const std::string &output_filename, double alpha, double to_the_power_of,
    std::string subversion, Logger &logger) {
    std::vector<Segment_w_info> input_segments;

    read_in(input_segments, input_filename, logger);

    logger.start_operation();
    Grid grid = SUBDIVISION::root_grid(input_segments, to_the_power_of);
    std::vector<std::vector<Segment_w_info>> subdivision;
    if (subversion=="0") {
        subdivision = SUBDIVISION::subdivide_grid(input_segments, grid);
    }
    else if (subversion=="1") {
        subdivision = SUBDIVISION::subdivide_grid_recursive(input_segments, grid, to_the_power_of);
    }
    else{
        throw std::runtime_error("subversion not recognized");
    }
    SUBDIVISION::plot_grid(input_segments,grid, "grid.svg");


    logger.end_operation("Subdivision (milliseconds)");
    logger.add("Number of subdivisions", subdivision.size());

    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(input_segments),"input.svg");

    int64_t extending_time = 0, build_arr_time = 0, build_graph_time = 0, max_flow_time = 0, combining_time = 0;
    std::vector<Polygon_2> polygons(0);

    for (auto &subset : subdivision) {

        if (subset.empty()) continue;
        EDGE_EXTENSION::add_outer_box(subset, 0.01);

        logger.start_operation();
        std::vector<Segment_w_info> extended_segments = EDGE_EXTENSION::STANDARD::extension(subset);
        extending_time += logger.operation_duration().count();


        logger.start_operation();
        Arrangement arr = ARRANGEMENT::build_arrangement(extended_segments, logger);
        build_arr_time += logger.operation_duration().count();


        logger.start_operation();
        Graph graph = GRAPH::build_graph(arr, alpha, logger);
        build_graph_time += logger.operation_duration().count();


        logger.start_operation();
        std::vector<bool> max_flow_solution = MAX_FLOW::max_flow(graph, arr);
        max_flow_time += logger.operation_duration().count();

        for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
            if (max_flow_solution[fit->data().id]) {
                Polygon_2 poly;
                auto curr = fit->outer_ccb();
                do {
                    poly.push_back(curr->source()->point());
                    ++curr;
                } while (curr != fit->outer_ccb());
                polygons.emplace_back(poly);
            }
        }
    }
    logger.add("Subdivided extending (milliseconds)", extending_time);
    logger.add("Subdivided build arr (milliseconds)", build_arr_time);
    logger.add("Subdivided build graph (milliseconds)", build_graph_time);
    logger.add("Subdivided max_flow (milliseconds)", max_flow_time);


    logger.start_operation();
    auto output_data = IO_FUNCTIONS::cgal_combines(polygons);
    logger.end_operation("Combining faces of solution (milliseconds)");

    IO_FUNCTIONS::GPKG::write_to_gpkg(output_data, logger.out_dir_stem() + "_subdivision");

    std::vector <Segment_w_info> segs;
    IO_FUNCTIONS::pwh_to_swi(output_data, segs);

    EDGE_EXTENSION::add_outer_box(segs, 0.01);

    logger.start_operation();
    auto extended = EDGE_EXTENSION::STANDARD::extension(segs);
    logger.end_operation("Final extending edges (milliseconds) ");

    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(extended),"extension_after_merge.svg");

    logger.start_operation();
    Arrangement arr = ARRANGEMENT::build_arrangement(extended, logger);
    logger.end_operation("Final build arr (milliseconds) ");

    logger.start_operation();
    Graph graph = GRAPH::build_graph(arr, alpha, logger);
    logger.end_operation("Final build graph (milliseconds) ");

    logger.start_operation();

    std::vector<bool> max_flow_solution = MAX_FLOW::max_flow(graph, arr);
    logger.end_operation("Final max flow (milliseconds) ");

    logger.start_operation();
    auto holes_and_outer = IO_FUNCTIONS::combine_polygons(max_flow_solution, arr);
    auto final_output = IO_FUNCTIONS::create_polygons_with_holes(holes_and_outer.first, holes_and_outer.second);
    logger.end_operation("Final combining faces of solution (milliseconds)");

    IO_FUNCTIONS::GPKG::write_to_gpkg(final_output, logger.out_dir_stem() + "_solution");

    std::string command = "qgis " + logger.out_dir_stem() + "_solution.gpkg " + logger.out_dir_stem() + "_subdivision.gpkg "
    + input_filename + " &";

    int status = std::system(command.c_str());


    logger.end();
}


void run_preprocessed(const std::string &input_filename, const std::string &output_filename, double alpha, double degree,
    double distance, std::string subversion, Logger &logger) {

    std::vector<Segment_w_info> input_segments;

    read_in(input_segments, input_filename, logger);

    EDGE_EXTENSION::add_outer_box(input_segments, 0.01);

    logger.start_operation();
    std::vector<Segment_w_info> extended_segments = EDGE_EXTENSION::STANDARD::extension(input_segments);
    logger.end_operation("Extending edges (milliseconds) ");


    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(extended_segments), "after_extension.svg");


    std::vector<std::vector<Segment_w_info>> spatially_close;
    if (false) {
        auto segs = PRE_PROCESS::group_degree(extended_segments, degree);

        spatially_close = PRE_PROCESS::spatially_close_groups(segs,
            distance);
    }
    else {
        spatially_close = PRE_PROCESS::group_by_degree_and_closeness(extended_segments, degree,distance);
    }
    if (subversion=="0") {
        PRE_PROCESS::longest_wins(spatially_close);
    }
    else if (subversion=="1") {
        PRE_PROCESS::shortest_wins(spatially_close);
    }
    else if (subversion=="2") {
        PRE_PROCESS::longest_and_shortest_wins(spatially_close);
    }
    else if (subversion=="3") {
        PRE_PROCESS::longest_mid_shortest_wins(spatially_close);
    }

    assert(spatially_close.size()>0);
    std::vector<Segment_w_info> filtered_extended = spatially_close[0];

    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(filtered_extended), "after_filtering.svg");

    std::vector<PWH> output_data(0);

    logger.start_operation();
    Arrangement arr = ARRANGEMENT::build_arrangement(extended_segments, logger);
    logger.end_operation("Building Arragnement (milliseconds)");
    logger.add("Number of edges in arrangement", arr.number_of_edges());
    logger.add("Number of halfedges in arrangement", arr.number_of_halfedges());
    logger.add("Number of faces in arrangement", arr.number_of_faces());


    aggregate(output_data, arr, logger);

    IO_FUNCTIONS::GPKG::write_to_gpkg(output_data, logger.out_dir_stem() + "_solution");

    std::string command = "qgis " + logger.out_dir_stem() + "_solution.gpkg " + input_filename + " &";
    int status = std::system(command.c_str());

    logger.end();
}



}
