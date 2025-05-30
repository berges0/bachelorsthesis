#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <limits>

#include "edge_extension.hpp"
#include "instance_loader/BasicDataStructures.hpp"
#include "instance_loader/Geometry.hpp"
#include "instance_loader/shapefile_io_operations.h"


void exportToSVG(const std::vector<SHPLoader::Point2D>& points, const std::vector<std::vector<int>>& polygons, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Could not write to " << filename << std::endl;
        return;
    }

    double width = 1000;
    double height = 1000;

    // Find bounding box of points
    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();

    for (const auto& p : points) {
        if (p.x < min_x) min_x = p.x;
        if (p.x > max_x) max_x = p.x;
        if (p.y < min_y) min_y = p.y;
        if (p.y > max_y) max_y = p.y;
    }

    // Function to normalize coordinate to SVG space
    auto normalize_x = [&](double x) {
        return (x - min_x) / (max_x - min_x) * width;
    };

    auto normalize_y = [&](double y) {
        // SVG y=0 is top, so invert y-axis
        return height - ((y - min_y) / (max_y - min_y) * height);
    };


    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << " " << height << "\">" << std::endl;

    for (const auto& poly : polygons) {
        out << "<polygon points=\"";
        for (int idx : poly) {
            if (idx < 0 || idx >= points.size()) continue;
            out << normalize_x(points[idx].x) << "," << normalize_y(points[idx].y) << " ";
        }
        out << "\" style=\"fill:none;stroke:black;stroke-width:1\" />" << std::endl;
    }

    out << "</svg>" << std::endl;
    out.close();
}

int main(int argc, char* argv[]) {
    std::cout << "Hello from Edge-Aligned Aggregation!" << std::endl;

    // Beispielcode vom Original-Repo einbinden
    auto pa = SHPLoader::Point2D(3.,4.);
    auto pb = SHPLoader::Point2D(3.,5.);
    std::cout << SHPLoader::distance(pa,pb) << std::endl;

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
        exportToSVG(data.first, data.second, "instance.svg");
        EdgeExtension Extension;
        Extension.buildRTree(data.first, data.second);

        auto nearest = Extension.queryNearest(0.5, 0.0, 3);
        for (const auto& seg : nearest) {
            std::cout << "Segment: " << bg::wkt<Segment>(seg) << "\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Exception while loading instance: " << e.what() << std::endl;
        return 3;
    }
    return 0;
}
