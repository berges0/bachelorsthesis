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
           version_input, clParser.getAlpha());

    if (version == "1" ) {
        Logger logger_1(clParser.inputFileName(), clParser.outputDirectory(),clParser.getTimeLimit().count(),
        version_input,clParser.getThreshold(), clParser.getThresholdScale(), clParser.getThresholdVariant());
        logger=logger_1;
    }

    logger.add("Input file name", clParser.inputFileName());
    logger.add("Output directory", clParser.outputDirectory());
    logger.add("Version", version_input);
    logger.add("Threshold", clParser.getThreshold());
    logger.add("Threshold scale", clParser.getThresholdScale());
    logger.add("Threshold_deg", clParser.getDegree());
    logger.add("Threshold_dist", clParser.getDistance());
    logger.add("To the power of", clParser.getPower());
    logger.add("Timelimit", clParser.getTimeLimit().count());
    logger.add("Threshold variant", clParser.getThresholdVariant());

    /*std::string shp_input;
    std::string json_string;
    std::string gpkg_ouput;
    for (const auto& entry : std::filesystem::directory_iterator(clParser.inputFileName())) {
        if (entry.is_regular_file()) {
            std::filesystem::path file_path = entry.path();
            std::string ext = file_path.extension().string();

            if (ext == ".shp") {
                shp_input = file_path.parent_path().string()+"/"+file_path.stem().string();
                gpkg_ouput= file_path.parent_path().string()+"/"+file_path.stem().string()+".gpkg";
            }
            else if (ext == ".json") {
                json_string = file_path.string();
            }
        }
    }

    std::vector<PWH> polys = IO_FUNCTIONS::SHP::read_shp_to_pwh(shp_input);

    // --- Compute area & perimeter ---
    Kernel::FT area(0);
    Kernel::FT perimeter(0);

    for (auto &pwh : polys) {
        // outer boundary
        area += CGAL::abs(pwh.outer_boundary().area());
        for (auto eit = pwh.outer_boundary().edges_begin(); eit != pwh.outer_boundary().edges_end(); ++eit) {
            perimeter += std::sqrt(CGAL::to_double(eit->squared_length()));
        }
        // holes
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
            std::cout<<"HAS HOLE";
            area -= CGAL::abs(hit->area());
            for (auto eit = hit->edges_begin(); eit != hit->edges_end(); ++eit) {
                perimeter += std::sqrt(CGAL::to_double(eit->squared_length()));
            }
        }
    }

    double d_area = CGAL::to_double(area);
    double d_perimeter = CGAL::to_double(perimeter);
    std::ifstream in(json_string);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open json file " + json_string);
    }

    nlohmann::json j;
    in >> j; // parse JSON
    in.close();

    double alpha = j["Alpha"];
    j["Area"]=d_area;
    j["Perimeter"]=d_perimeter;
    j["Objective Value"]=alpha*d_area + (1-alpha)*d_perimeter;
    // Modify or add an attribute

    // Write JSON back to file (overwrite)
    std::ofstream out(json_string);
    out << j.dump(4); // pretty print with indent of 4 spaces
    out.close();

    std::cout<< " PERIMETER: " << perimeter << '\n';
    std::cout<< " AREA: " << area << '\n';


    IO_FUNCTIONS::GPKG::write_to_gpkg(polys, gpkg_ouput);*/

    if (version == "0") {
        ALGORITHM::run_standard(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), logger);
    } else if (version == "1") {
        ALGORITHM::run_limited(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getThreshold(),
            clParser.getThresholdScale(), clParser.getThresholdVariant(), logger);
    } else if (version == "2"){
        ALGORITHM::run_preprocessed(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getDegree(),
            clParser.getDistance(), subversion, logger);
    } else if (version == "3") {
        ALGORITHM::run_edge_relink(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getThreshold(),
            clParser.getThresholdScale(), clParser.getThresholdVariant(), clParser.getDegree(), clParser.getDistance(), subversion, logger);
    } else if (version == "4") {
        ALGORITHM::run_subdivision(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getPower(),
            subversion, logger);
    } else if (version == "5") {
        ALGORITHM::run_outer_endpoints(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(),clParser.getThreshold(),
            clParser.getThresholdScale(), clParser.getThresholdVariant(), clParser.getDegree(), clParser.getDistance(), subversion, logger);
    } else{
        std::cerr << "Unknown version: " << version << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

