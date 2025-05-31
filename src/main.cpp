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

void write_segments_to_svg(const std::vector<Segment>& segments, const std::string& filename) {
    std::ofstream svg(filename);
    if (!svg.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // SVG header
    svg << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n";

    // Optionally, add a white background rectangle (for better visibility)
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\" />\n";

    // Draw each segment
    for (const auto& seg : segments) {
        double x1 = CGAL::to_double(seg.source().x());
        double y1 = CGAL::to_double(seg.source().y());
        double x2 = CGAL::to_double(seg.target().x());
        double y2 = CGAL::to_double(seg.target().y());

        // SVG Y-axis goes down, so you might want to flip Y if needed
        // Here, I assume coordinates are suitable as-is.
        svg << "<line x1=\"" << x1 << "\" y1=\"" << y1
            << "\" x2=\"" << x2 << "\" y2=\"" << y2
            << "\" stroke=\"black\" stroke-width=\"1\" />\n";
    }

    svg << "</svg>\n";
    svg.close();
    std::cout << "SVG written to " << filename << std::endl;
}


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
        write_segments_to_svg(algo.build_segment_vector(data.first, data.second),"original.svg");
        write_segments_to_svg(algo.ray_shoot_intersection(algo.build_segment_vector(data.first, data.second)),"extended.svg");


    } catch (const std::exception& e) {
        std::cerr << "Exception while loading instance: " << e.what() << std::endl;
        return 3;
    }
    return 0;
}
