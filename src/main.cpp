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

    Logger logger(clParser.jsonFileName());

    logger.add("Input file name", clParser.inputFileName());
    logger.add("Output file name", clParser.outputFileName());
    logger.add("Json Filename", clParser.jsonFileName());
    logger.add("Alpha", clParser.getAlpha());
    logger.add("Version", version_input);
    logger.add("Threshold", clParser.getThreshold());
    logger.add("Threshold_deg", clParser.getDegree());
    logger.add("Threshold_dist", clParser.getDistance());
    logger.add("To the power of", clParser.getPower());
    logger.add("Timelimit", clParser.getTimeLimit().count());

    std::string version = "0";
    std::string subversion = "0";

    if (version_input.size()>1) {
        assert(version_input.size()==3);
        version = version_input.substr(0);
        subversion = version_input.substr(3, version_input.size()-1);
    }
    else if (version_input.size()==1){
        version = version_input;
    }

    if (version == "0") {
        ALGORITHM::run_standard(clParser.inputFileName(), clParser.outputFileName(), clParser.getAlpha(), logger);
    } else if (version == "1") {
        ALGORITHM::run_limited(clParser.inputFileName(), clParser.outputFileName(), clParser.getAlpha(), clParser.getThreshold(),
            logger);
    } else if (version == "2"){
        ALGORITHM::run_preprocessed(clParser.inputFileName(), clParser.outputFileName(), clParser.getAlpha(), clParser.getDegree(),
            clParser.getDistance(), subversion, logger);
    } else if (version == "3") {
        ALGORITHM::run_edge_relink(clParser.inputFileName(), clParser.outputFileName(), clParser.getAlpha(), clParser.getDegree(),
            clParser.getDistance(), subversion, logger);
    } else if (version == "4") {
        ALGORITHM::run_subdivision(clParser.inputFileName(), clParser.outputFileName(), clParser.getAlpha(), clParser.getPower(),
            subversion, logger);
    }
    else{
        std::cerr << "Unknown version: " << version << ". Supported versions are: standard, limited, subdivision.\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

