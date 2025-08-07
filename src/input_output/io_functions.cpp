#include "algorithm.hpp"

void algorithm::segments_to_svg(const std::vector<Segment>& segments, const std::string& filename) {

    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    for (const auto& s : segments) {
        minX = std::min({minX, CGAL::to_double(s.source().x()), CGAL::to_double(s.target().x())});
        maxX = std::max({maxX, CGAL::to_double(s.source().x()), CGAL::to_double(s.target().x())});
        minY = std::min({minY, CGAL::to_double(s.source().y()), CGAL::to_double(s.target().y())});
        maxY = std::max({maxY, CGAL::to_double(s.source().y()), CGAL::to_double(s.target().y())});
    }

    double padding = 10.0; // Adjust based on your coordinate range
    minX -= padding;
    maxX += padding;
    minY -= padding;
    maxY += padding;

    double width = maxX - minX;
    double height = maxY - minY;

    std::ofstream svg_file(filename);
    svg_file << std::fixed << std::setprecision(10); // HIGH precision output
    double pixel_width = 1000.0;
    double aspect_ratio = height / width;
    double pixel_height = pixel_width * aspect_ratio;
    double scale = pixel_width / width; // pixel per unit
    double stroke_width_units = 1.0 / scale; // 1 screen pixel = X world units

    svg_file << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
             << "viewBox=\"" << minX << " " << minY << " " << width << " " << height << "\" "
             << "width=\"" << pixel_width << "\" height=\"" << pixel_height << "\" "
             << "preserveAspectRatio=\"xMidYMid meet\" "
             << "fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width_units << "\">\n";


    for (const auto& s : segments) {
        svg_file << "<line x1=\"" << s.source().x() << "\" y1=\"" << (maxY -(s.source().y()- minY))
                 << "\" x2=\"" << s.target().x() << "\" y2=\"" << (maxY - (s.target().y() - minY))
                 << "\" />\n";
    }

    svg_file << "</svg>\n";

    svg_file.close();
}
void algorithm::polygons_to_svg(const std::vector<Polygon_with_holes_2>& polygons, const std::string& filename) {
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    // Compute bounding box over all outer and hole polygons
    for (const auto& pwh : polygons) {
        for (const auto& p : pwh.outer_boundary()) {
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            minX = std::min(minX, x);
            maxX = std::max(maxX, x);
            minY = std::min(minY, y);
            maxY = std::max(maxY, y);
        }
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
            for (const auto& p : *hit) {
                double x = CGAL::to_double(p.x());
                double y = CGAL::to_double(p.y());
                minX = std::min(minX, x);
                maxX = std::max(maxX, x);
                minY = std::min(minY, y);
                maxY = std::max(maxY, y);
            }
        }
    }

    double padding = 10.0;
    minX -= padding;
    maxX += padding;
    minY -= padding;
    maxY += padding;

    double width = maxX - minX;
    double height = maxY - minY;

    std::ofstream svg_file(filename);
    svg_file << std::fixed << std::setprecision(10);

    double pixel_width = 1000.0;
    double aspect_ratio = height / width;
    double pixel_height = pixel_width * aspect_ratio;
    double scale = pixel_width / width;
    double stroke_width_units = 1.0 / scale;

    svg_file << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
             << "viewBox=\"" << minX << " " << minY << " " << width << " " << height << "\" "
             << "width=\"" << pixel_width << "\" height=\"" << pixel_height << "\" "
             << "preserveAspectRatio=\"xMidYMid meet\" "
             << "fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width_units << "\">\n";

    for (const auto& pwh : polygons) {
        // Outer boundary
        svg_file << "<polygon points=\"";
        for (const auto& p : pwh.outer_boundary()) {
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            svg_file << x << "," << (maxY - (y - minY)) << " ";
        }
        svg_file << "\" fill=\"lightgray\" stroke=\"black\"/>\n";

        // Holes
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
            svg_file << "<polygon points=\"";
            for (const auto& p : *hit) {
                double x = CGAL::to_double(p.x());
                double y = CGAL::to_double(p.y());
                svg_file << x << "," << (maxY - (y - minY)) << " ";
            }
            svg_file << "\" fill=\"white\" stroke=\"black\"/>\n";
        }
    }

    svg_file << "</svg>\n";
    svg_file.close();
}

void algorithm::polygonns_to_svg(const std::vector<Polygon_2>& polygons, const std::string& filename) {
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    // Compute bounding box over all polygons
    for (const auto& poly : polygons) {
        for (const auto& p : poly) {
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            minX = std::min(minX, x);
            maxX = std::max(maxX, x);
            minY = std::min(minY, y);
            maxY = std::max(maxY, y);
        }
    }

    double padding = 10.0;
    minX -= padding;
    maxX += padding;
    minY -= padding;
    maxY += padding;

    double width = maxX - minX;
    double height = maxY - minY;

    std::ofstream svg_file(filename);
    svg_file << std::fixed << std::setprecision(10);

    double pixel_width = 1000.0;
    double aspect_ratio = height / width;
    double pixel_height = pixel_width * aspect_ratio;
    double scale = pixel_width / width;
    double stroke_width_units = 1.0 / scale;

    svg_file << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
             << "viewBox=\"" << minX << " " << minY << " " << width << " " << height << "\" "
             << "width=\"" << pixel_width << "\" height=\"" << pixel_height << "\" "
             << "preserveAspectRatio=\"xMidYMid meet\" "
             << "fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width_units << "\">\n";

    // Draw each polygon
    for (const auto& poly : polygons) {
        svg_file << "<polygon points=\"";
        for (const auto& p : poly) {
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            // Flip Y-axis for SVG
            svg_file << x << "," << (maxY - (y - minY)) << " ";
        }
        svg_file << "\" />\n";
    }

    svg_file << "</svg>\n";
    svg_file.close();
}
void algorithm::write_to_shp(const std::vector<Polygon_with_holes_2>& polygons, const std::string& filename) {
    std::vector<std::pair<double, double>> points;
    std::vector<std::vector<int>> polys;

    std::vector<int> polygon_indices;

    for (const auto& pwh : polygons) {
        // Outer boundary
        polygon_indices.clear();
        for (const auto& p : pwh.outer_boundary()) {
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            points.emplace_back(x, y);
            polygon_indices.emplace_back(points.size() - 1);
        }
        polys.push_back(polygon_indices);

        // Holes
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
            polygon_indices.clear();
            for (const auto& p : *hit) {
                double x = CGAL::to_double(p.x());
                double y = CGAL::to_double(p.y());
                points.emplace_back(x, y);
                polygon_indices.emplace_back(points.size() - 1);
            }
            polys.push_back(polygon_indices);
        }
    }

    std::cout << "THERE ARE " << points.size() << " points" << std::endl;
    std::cout << "AND " << polys.size() << " polygons (including holes)" << std::endl;

    SHPLoader::writeToShapeFile(std::make_pair(points, polys), filename);
}
