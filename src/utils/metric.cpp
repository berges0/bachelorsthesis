#include "utils/metric.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "io/io_functions.hpp"

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: metric <INPUT.txt> <OUTPUT.txt> <METRIC_VERSION>\n";
        return 1;
    }

    const std::string input  = argv[1];
    const std::string output = argv[2];
    const std::string version = argv[3];

    std::vector<PWH> inputpolys = IO_FUNCTIONS::GPKG::read_gpkg_to_pwh(input);
    std::vector<PWH> outputpolys = IO_FUNCTIONS::GPKG::read_gpkg_to_pwh(output);

    if (version )




}
