//
// Created by samuel-berge on 8/11/25.
//
#include <fstream>
#include <iostream>
#include <vector>
#include "io/chrono_io.hpp"
#include "io/cl_parser.hpp"
#include "io/io_functions.hpp"
#include "core/algorithm.hpp"
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

    Logger logger("output_1.json");

    logger.add("Input file name", clParser.inputFileName());
    logger.add("Output file name", clParser.outputFileName());
    logger.add("Alpha", alpha);
    logger.add("Version", version);
    logger.add("Threshold", threshold);
    logger.add("Timelimit", clParser.getTimeLimit().count());

    if (version == "0") {
        ALGORITHM::run_standard(clParser.inputFileName(), clParser.outputFileName(), alpha, logger);
    } else if (version == "1") {
        ALGORITHM::run_limited(clParser.inputFileName(), clParser.outputFileName(), alpha, threshold, logger);
    } else if (version == "2") {
        ALGORITHM::run_subdivision(clParser.inputFileName(), clParser.outputFileName(), alpha, threshold, logger);
    } else {
        std::cerr << "Unknown version: " << version << ". Supported versions are: standard, limited, subdivision.\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

