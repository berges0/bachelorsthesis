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

    std::string version_input = clParser.getVersion();

    std::string version = "0";
    std::string subversion = "0";

    if (version_input.size()>1) {
        assert(version_input.size()==3);
        version = version_input.substr(0,1);
        subversion = version_input.substr(2, 1);
    }
    else if (version_input.size()==1){
        version = version_input;
    }

    Logger logger(clParser.inputFileName(), clParser.outputDirectory(),clParser.getTimeLimit().count(),
        version_input);

    logger.add("Input file name", clParser.inputFileName());
    logger.add("Output directory", clParser.outputDirectory());
    logger.add("Alpha", clParser.getAlpha());
    logger.add("Version", version_input);
    logger.add("Threshold scale", clParser.getThresholdScale());
    logger.add("Threshold_deg", clParser.getDegree());
    logger.add("Threshold_dist", clParser.getDistance());
    logger.add("To the power of", clParser.getPower());
    logger.add("Timelimit", clParser.getTimeLimit().count());
    logger.add("Threshold variant", clParser.getThresholdVariant());


    if (version == "0") {
        ALGORITHM::run_standard(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), logger);
    } else if (version == "1") {
        ALGORITHM::run_limited(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getThresholdScale(),
            clParser.getThresholdVariant(), logger);
    } else if (version == "2"){
        ALGORITHM::run_preprocessed(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getDegree(),
            clParser.getDistance(), subversion, logger);
    } else if (version == "3") {
        ALGORITHM::run_edge_relink(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getThresholdScale(),
        clParser.getThresholdVariant(), clParser.getDegree(), clParser.getDistance(), subversion, logger);
    } else if (version == "4") {
        ALGORITHM::run_subdivision(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getPower(),
            subversion, logger);
    } else if (version == "5") {
        ALGORITHM::run_outer_endpoints(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(),clParser.getThresholdScale(),
            clParser.getThresholdVariant(), clParser.getDegree(), clParser.getDistance(), subversion, logger);
    } else{
        std::cerr << "Unknown version: " << version << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

