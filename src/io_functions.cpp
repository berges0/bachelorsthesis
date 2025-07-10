#include "algorithm.hpp"

void algorithm::segments_to_svg(const std::vector<Segment>& segments, const std::string& filename) {

    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();

    for (const auto& s : segments) {
        min_x = std::min({min_x, CGAL::to_double(s.source().x()), CGAL::to_double(s.target().x())});
        max_x = std::max({max_x, CGAL::to_double(s.source().x()), CGAL::to_double(s.target().x())});
        min_y = std::min({min_y, CGAL::to_double(s.source().y()), CGAL::to_double(s.target().y())});
        max_y = std::max({max_y, CGAL::to_double(s.source().y()), CGAL::to_double(s.target().y())});
    }

    double padding = 10.0; // Adjust based on your coordinate range
    min_x -= padding;
    max_x += padding;
    min_y -= padding;
    max_y += padding;

    double width = max_x - min_x;
    double height = max_y - min_y;

    std::ofstream svg_file(filename);
    svg_file << std::fixed << std::setprecision(10); // HIGH precision output
    double pixel_width = 1000.0;
    double aspect_ratio = height / width;
    double pixel_height = pixel_width * aspect_ratio;
    double scale = pixel_width / width; // pixel per unit
    double stroke_width_units = 1.0 / scale; // 1 screen pixel = X world units

    svg_file << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
             << "viewBox=\"" << min_x << " " << min_y << " " << width << " " << height << "\" "
             << "width=\"" << pixel_width << "\" height=\"" << pixel_height << "\" "
             << "preserveAspectRatio=\"xMidYMid meet\" "
             << "fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width_units << "\">\n";


    for (const auto& s : segments) {
        svg_file << "<line x1=\"" << s.source().x() << "\" y1=\"" << (max_y -(s.source().y()- min_y))
                 << "\" x2=\"" << s.target().x() << "\" y2=\"" << (max_y - (s.target().y() - min_y))
                 << "\" />\n";
    }

    svg_file << "</svg>\n";

    svg_file.close();
}
