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

/*
void algorithm::poly_to_svg(const std::vector<Polygon_2>& polygons, const std::string& filename) {

    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    for (const auto& poly : polygons) {
        for (const Point& p : poly) {
            if (p.x() < minX) minX = CGAL::to_double(p.x());
            if (p.x() > maxX) maxX = CGAL::to_double(p.x());
            if (p.y() < minY) minY = CGAL::to_double(p.y());
            if (p.y() > maxY) maxY = CGAL::to_double(p.y());
        }
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



*/