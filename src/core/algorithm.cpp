//
// Created by samuel-berge on 8/20/25.
//
#include "core/algorithm.hpp"

#include "core/edge_relink.hpp"

namespace ALGORITHM{

void read_in (std::vector<Segment_w_info> &input_segments, const std::string &input_filename, Logger &logger) {
    auto data = SHPLoader::ReadShapeFileToPoint2D(input_filename);
    logger.add("Number of input polygons", data.second.size());

    input_segments= IO_FUNCTIONS::read_in_segment(data.first, data.second);
    logger.add("Number of input segments", input_segments.size());

}

void aggregate(std::vector<Polygon_with_holes_2> &output_data, const std::vector<Segment_w_info> &input, double alpha, Logger &logger) {

    logger.start_operation();
    Arrangement arr = ARRANGEMENT::build_arrangement(input);
    logger.end_operation("Building Arragnement (milliseconds)");
    logger.add("Number of segments after extension", arr.number_of_edges());
    logger.add("Number of faces in arrangement", arr.number_of_faces());

    logger.start_operation();
    Graph graph = GRAPH::build_graph(arr, alpha);
    logger.end_operation("Building Graph (milliseconds)");
    logger.add("Number of edges in graph", std::get<0>(graph).size());

    logger.start_operation();
    std::vector<bool> max_flow_solution = MAX_FLOW::max_flow(graph, arr);
    logger.end_operation("Running Max Flow (milliseconds)");
    int nr_faces_solution = 0;
    for (int i = 0; i < max_flow_solution.size(); i++) {if (max_flow_solution[i]) nr_faces_solution++;}
    logger.add("Number of faces solution", nr_faces_solution);

    logger.start_operation();
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

    //JUST FOR DEBUGGING
    IO_FUNCTIONS::segments_to_svg(EDGE_EXTENSION::filter_segments(extended_segments),"extended.svg");

    std::vector<Polygon_with_holes_2> output_data;

    aggregate(output_data, extended_segments, alpha, logger);

    IO_FUNCTIONS::writeToShapeFile(output_data, output_filename);

    std::string command1 = "qgis " + input_filename + " &";
    std::string command2 = "qgis " + output_filename + " &";
    int i = std::system(command1.c_str());
    int j = std::system(command2.c_str());
    logger.end();
}

void run_limited(const std::string &input_filename, const std::string &output_filename, double alpha, double threshold,
    Logger &logger) {

    std::vector<Segment_w_info> input_segments;

    read_in(input_segments, input_filename, logger);

    EDGE_EXTENSION::add_outer_box(input_segments,0.01);

    logger.start_operation();
    std::vector<Segment_w_info> extended_segments = EDGE_EXTENSION::LIMITED::extension(input_segments, threshold);
    logger.end_operation("Extending edges (milliseconds) ");

    std::vector<Polygon_with_holes_2> output_data;

    aggregate(output_data, extended_segments, alpha, logger);

    IO_FUNCTIONS::writeToShapeFile(output_data, output_filename);

    std::string command1 = "qgis " + input_filename + " &";
    std::string command2 = "qgis " + output_filename + " &";
    int i = std::system(command1.c_str());
    int j = std::system(command2.c_str());
    logger.end();
}

void run_preprocessed(const std::string &input_filename, const std::string &output_filename, double alpha, double degree,
    double distance, std::string subversion, Logger &logger) {

    auto data = SHPLoader::ReadShapeFileToPoint2D(input_filename);
    logger.add("Number of input polygons", data.second.size());

    std::vector<Segment_w_info> input_segments= IO_FUNCTIONS::read_in_segment(data.first, data.second);
    logger.add("Number of input segments", input_segments.size());

    EDGE_EXTENSION::add_outer_box(input_segments, 0.01);

    logger.start_operation();
    std::vector<Segment_w_info> extended_segments = EDGE_EXTENSION::STANDARD::extension(input_segments);
    logger.end_operation("Extending edges (milliseconds) ");

    IO_FUNCTIONS::segments_to_svg(EDGE_EXTENSION::filter_segments(extended_segments), "after_extension.svg");

    auto segs = PRE_PROCESS::group_degree(extended_segments, degree);

    auto spatially_close = PRE_PROCESS::spatially_close_groups(segs,
        distance);

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

    std::vector<Polygon_with_holes_2> output_data(0);

    aggregate(output_data, filtered_extended, alpha, logger);

    IO_FUNCTIONS::writeToShapeFile(output_data, output_filename);

    std::string command1 = "qgis " + input_filename + " &";
    std::string command2 = "qgis " + output_filename + " &";
    int i = std::system(command1.c_str());
    int j = std::system(command2.c_str());
    logger.end();
}

void run_edge_relink(const std::string &input_filename, const std::string &output_filename, double alpha, double degree,
                     double distance, std::string subversion, Logger &logger) {

    std::vector<Segment_w_info> input_segments;

    read_in(input_segments, input_filename, logger);

    EDGE_EXTENSION::add_outer_box(input_segments,0.01);

    logger.start_operation();
    std::vector<Segment_w_info> extended_segments = EDGE_EXTENSION::STANDARD::extension(input_segments);
    logger.end_operation("Extending edges (milliseconds) ");

    auto segs = PRE_PROCESS::group_degree(extended_segments, degree);

    auto spatially_close = PRE_PROCESS::spatially_close_groups(segs,
        distance);

    EDGE_RELINK::relink_edges(spatially_close);

    auto relinked_segments = spatially_close[0];

    std::vector<Polygon_with_holes_2> output_data(0);

    aggregate(output_data, relinked_segments, alpha, logger);

    IO_FUNCTIONS::writeToShapeFile(output_data, output_filename);

    std::string command1 = "qgis " + input_filename + " &";
    std::string command2 = "qgis " + output_filename + " &";
    int i = std::system(command1.c_str());
    int j = std::system(command2.c_str());
    logger.end();
}


void run_subdivision(const std::string &input_filename, const std::string &output_filename, double alpha, double to_the_power_of,
    std::string subversion, Logger &logger) {
    auto data = SHPLoader::ReadShapeFileToPoint2D(input_filename);
    logger.add("Number of input polygons", data.second.size());

    auto input_segments = IO_FUNCTIONS::read_in_segment(data.first, data.second);
    logger.add("Number of input segments", input_segments.size());

    logger.start_operation();
    Grid grid = SUBDIVISION::root_grid(input_segments, 0.77);
    auto subdivision = SUBDIVISION::subdivide_grid(input_segments,data.second.size(), grid);
    logger.end_operation("Subdivision (milliseconds)");
    logger.add("Number of subdivisions", subdivision.size());

    IO_FUNCTIONS::segments_to_svg(EDGE_EXTENSION::filter_segments(input_segments),"input.svg");
    SUBDIVISION::plot_grid(input_segments, grid);

    int64_t extending_time = 0, build_arr_time = 0, build_graph_time = 0, max_flow_time = 0, combining_time = 0;
    std::vector<Polygon_2> polygons(0);

    for (auto &subset : subdivision) {
        if (subset.empty()) continue;
        EDGE_EXTENSION::add_outer_box(subset, 0.01);

        logger.start_operation();
        std::vector<Segment_w_info> extended_segments = EDGE_EXTENSION::STANDARD::extension(subset);
        extending_time += logger.operation_duration().count();

        logger.start_operation();
        Arrangement arr = ARRANGEMENT::build_arrangement(extended_segments);
        build_arr_time += logger.operation_duration().count();

        logger.start_operation();
        Graph graph = GRAPH::build_graph(arr, alpha);
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
    assert (output_filename.size() >= 4 && output_filename.compare(output_filename.size() - 4, 4, ".shp") == 0);
    std::string subdiv_filename = output_filename.substr(0, output_filename.size() - 4) + "_subdivision.shp";

    IO_FUNCTIONS::writeToShapeFile(output_data,subdiv_filename);
    std::string command1 = "qgis " + subdiv_filename + " &";
    int i = std::system(command1.c_str());

    auto segs = SUBDIVISION::pwh_to_swi(output_data);

    EDGE_EXTENSION::add_outer_box(segs, 0.01);

    logger.start_operation();
    auto extended = EDGE_EXTENSION::STANDARD::extension(segs);
    logger.end_operation("Final extending edges (milliseconds) ");

    IO_FUNCTIONS::segments_to_svg(EDGE_EXTENSION::filter_segments(extended),"extension_after_merge.svg");

    logger.start_operation();
    Arrangement arr = ARRANGEMENT::build_arrangement(extended);
    logger.end_operation("Final build arr (milliseconds) ");

    logger.start_operation();
    Graph graph = GRAPH::build_graph(arr, alpha);
    logger.end_operation("Final build graph (milliseconds) ");

    logger.start_operation();
    std::vector<bool> max_flow_solution = MAX_FLOW::max_flow(graph, arr);
    logger.end_operation("Final max flow (milliseconds) ");

    logger.start_operation();
    auto holes_and_outer = IO_FUNCTIONS::combine_polygons(max_flow_solution, arr);
    auto final_output = IO_FUNCTIONS::create_polygons_with_holes(holes_and_outer.first, holes_and_outer.second);
    logger.end_operation("Final combining faces of solution (milliseconds)");

    IO_FUNCTIONS::writeToShapeFile(final_output, output_filename);
    std::string command2 = "qgis " + input_filename + " &";
    std::string command3 = "qgis " + output_filename + " &";
    int j = std::system(command2.c_str());
    int k = std::system(command3.c_str());

    logger.end();
}

}
