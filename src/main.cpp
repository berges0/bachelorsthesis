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

    Timer t;
    t.start();
    std::cout << "Running algorithm on input file name: " << clParser.inputFileName() << '\n';
    std::cout << "Alpha is: " << clParser.getAlpha() << '\n';
    std::cout << "Timelimit is: " << clParser.getTimeLimit() << '\n';
    std::cout << "Algorithm version is: " << clParser.getVersion() << '\n';

    std::string version = clParser.getVersion();
    double alpha = clParser.getAlpha();
    double threshold = 4; // Default threshold for limited version

    auto data = SHPLoader::ReadShapeFileToPoint2D(clParser.inputFileName());

    std::vector<Segment_w_info> input_segments= IO_FUNCTIONS::read_in_segment(data.first, data.second);

    EDGE_EXTENSION::add_outer_box(input_segments);

    std::vector<Segment_w_info> extended_segments = EDGE_EXTENSION::edge_extension(input_segments, version, threshold);

    IO_FUNCTIONS::segments_to_svg(EDGE_EXTENSION::filter_segments(extended_segments),"extended.svg");

    Arrangement arr = ARRANGEMENT::build_arrangement(extended_segments);

    Graph graph = GRAPH::build_graph(arr, alpha);

    std::vector<bool> max_flow_solution = MAX_FLOW::max_flow(graph, arr);

    std::vector<Polygon_with_holes_2> output_data = IO_FUNCTIONS::combine_polygons(max_flow_solution, arr);

    IO_FUNCTIONS::write_to_shp(output_data, clParser.outputFileName());

    std::string command1 = "qgis " + clParser.inputFileName() + " &";
    std::string command2 = "qgis " + clParser.outputFileName() + " &";
    int i = std::system(command1.c_str());
    int j = std::system(command2.c_str());

    std::cout << "Time taken overall: " << t.elapsed() << '\n';
    return EXIT_SUCCESS;
}

