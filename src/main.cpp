//
// Created by samuel-berge on 8/11/25.
//
#include <fstream>
#include <iostream>
#include <vector>
#include "io/chrono_io.hpp"
#include "io/cl_parser.hpp"
#include "io/io_functions.hpp"
#include "core/arrangement.hpp"
#include "core/edge_extension.hpp"
#include "core/graph.hpp"
#include "core/max_flow.hpp"
#include "utils/logger.hpp"

void printCommandLineFlags(int argc, char **argv, std::ostream &out = std::cout) {
    out << "Num args: " << argc << ".\nArgs: ";
    for (int i = 0; i < argc; i++) {
        out << argv[i] << ", ";
    }
    out << std::endl;
}

int main(int argc, char **argv) {
    printCommandLineFlags(argc, argv);
    CLParser clParser(argc, argv);

    std::cout << "Running algorithm on input file name: " << clParser.inputFileName() << '\n';
    std::cout << "Alpha is: " << clParser.getAlpha() << '\n';
    std::cout << "Timelimit is: " << clParser.getTimeLimit() << '\n';
    std::cout << "Algorithm version is: " << clParser.getVersion() << '\n';

    std::string version = clParser.getVersion();
    double alpha = clParser.getAlpha();
    double threshold = 4; // Default threshold for limited version

    Logger logger("output.json");

    logger.add("Input file name", clParser.inputFileName());
    logger.add("Output file name", clParser.outputFileName());
    logger.add("Alpha", alpha);
    logger.add("Version", version);
    logger.add("Threshold", threshold);
    logger.add("Timelimit", clParser.getTimeLimit().count());


    auto data = SHPLoader::ReadShapeFileToPoint2D(clParser.inputFileName());
    logger.add("Number of input polygons", data.second.size());

    std::vector<Segment_w_info> input_segments= IO_FUNCTIONS::read_in_segment(data.first, data.second);
    logger.add("Number of input segments", input_segments.size());

    EDGE_EXTENSION::add_outer_box(input_segments);

    logger.start_operation();
    std::vector<Segment_w_info> extended_segments = EDGE_EXTENSION::edge_extension(input_segments, version, threshold);
    logger.end_operation();
    logger.add("Extending edges (milliseconds)", logger.operation_duration().count());

    IO_FUNCTIONS::segments_to_svg(EDGE_EXTENSION::filter_segments(extended_segments),"extended.svg");

    logger.start_operation();
    Arrangement arr = ARRANGEMENT::build_arrangement(extended_segments);
    logger.end_operation();
    logger.add("Building Arragnement (milliseconds)", logger.operation_duration().count());
    logger.add("Number of segments after extension", arr.number_of_edges());
    logger.add("Number of faces in arrangement", arr.number_of_faces());

    logger.start_operation();
    Graph graph = GRAPH::build_graph(arr, alpha);
    logger.end_operation();
    logger.add("Building Graph (milliseconds)", logger.operation_duration().count());
    logger.add("Number of edges in graph", std::get<0>(graph).size());

    logger.start_operation();
    std::vector<bool> max_flow_solution = MAX_FLOW::max_flow(graph, arr);
    logger.end_operation();
    logger.add("Running Max Flow (milliseconds)", logger.operation_duration().count());
    int nr_faces_solution = 0;
    for (int i = 0; i < max_flow_solution.size(); i++) {if (max_flow_solution[i]) nr_faces_solution++;}
    logger.add("Number of faces solution", nr_faces_solution);

    logger.start_operation();
    auto holes_and_outer = IO_FUNCTIONS::combine_polygons(max_flow_solution, arr);
    auto output_data = IO_FUNCTIONS::create_polygons_with_holes(holes_and_outer.first, holes_and_outer.second);
    logger.end_operation();
    logger.add("Combining faces of solution (milliseconds)", logger.operation_duration().count());
    logger.add("Number of polygons in solution", holes_and_outer.first.size());
    logger.add("Number of holes in solution", holes_and_outer.second.size());

    IO_FUNCTIONS::writeToShapeFile(output_data, clParser.outputFileName());

    std::string command1 = "qgis " + clParser.inputFileName() + " &";
    std::string command2 = "qgis " + clParser.outputFileName() + " &";
    int i = std::system(command1.c_str());
    int j = std::system(command2.c_str());
    logger.end();

    return EXIT_SUCCESS;
}

