#include <algorithm.hpp>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "instance_loader/BasicDataStructures.hpp"
#include "instance_loader/Geometry.hpp"
#include "instance_loader/shapefile_io_operations.h"


#include <fstream>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;
typedef Kernel::Segment_2 Segment;


int main(int argc, char* argv[]) {
    std::cout << "Hello from Edge-Aligned Aggregation!" << std::endl;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_shapefile_path>" << std::endl;
        return 1;
    }

    std::string input_path = argv[1];

    // Check if file exists
    if (!std::filesystem::exists(input_path)) {
        std::cerr << "Error: Input file does not exist: " << input_path << std::endl;
        return 2;
    }

    try {
        // Load instance using the library
        auto data = SHPLoader::ReadShapeFileToPoint2D(input_path);
        std::cout << "Successfully loaded instance from: " << input_path << std::endl;
        std::cout << "Points: " << data.first.size() << std::endl;
        std::cout << "Polygons: " << data.second.size() << std::endl;
        algorithm algo;
        algo.segments_to_svg(algo.build_segment_vector(data.first, data.second),"original.svg");
        algo.segments_to_svg(algo.ray_shoot_intersection(algo.build_segment_vector(data.first, data.second)),"extended.svg");


    } catch (const std::exception& e) {
        std::cerr << "Exception while loading instance: " << e.what() << std::endl;
        return 3;
    }
    return 0;
}
