//
// Created by samuel-berge on 8/20/25.
//

#include "core/subdivision.hpp"

#include "core/edge_extension.hpp"
#include "io/io_functions.hpp"

namespace SUBDIVISION {

std::vector<std::vector<Segment_w_info>> subdivide_random(std::vector<Segment_w_info> &segments, double subset_size) {
    int n = segments.size();
    int k = log(n) * static_cast<int>(sqrt(n)); // chunk size = sqrt(n)
    std::vector<std::vector<Segment_w_info>> result;
    for (int i = 1; i < 3; ++i) {
        int index = 0;
        int limit = k / i;
        int subdivision_count = 1;
        while (limit < n) {
            std::vector<Segment_w_info> chunk;
            while (index < limit) {
                chunk.emplace_back(segments[index++]);
            }
            int lastpolyid = segments[index].poly_id;
            while (segments[index].poly_id == lastpolyid && index < n) {
                chunk.emplace_back(segments[index++]);
            }
            result.push_back(chunk);
            limit += k;
        }
        while (index < n) {
            std::vector<Segment_w_info> chunk;
            while (index < n) {
                chunk.emplace_back(segments[index++]);
            }
            result.push_back(chunk);
        }
    }
    return result;
}
Grid root_grid(const std::vector<Segment_w_info> &segments, double to_the_power_of) {
    Grid grid(segments);

    // nr of segments
    double n = static_cast<double>(grid.nr_segs);
    std::cout << "Number of segments: " << n << std::endl;

    // size of a subset should be s
    double s = std::pow(n, to_the_power_of);
    std::cout << "Size of a subset: " << s << std::endl;

    // number of subsets
    double m = std::ceil(n / s);
    std::cout << "Number of subsets: " << m << std::endl;

    // height/length ratio
    double a = grid.lenY / grid.lenX;

    grid.nr_cells_y = std::max(1., std::floor(std::sqrt(m * a)));
    std::cout << "Grid nr of cells in y: " << grid.nr_cells_y << std::endl;

    grid.nr_cells_x = std::ceil(m / grid.nr_cells_y);
    std::cout << "Grid nr of cells in x: " << grid.nr_cells_x << std::endl;

    grid.nr_cells = grid.nr_cells_x * grid.nr_cells_y;
    std::cout << "Grid nr of cells: " << grid.nr_cells << std::endl;

    grid.stepX = std::ceil(grid.lenX / grid.nr_cells_x);
    grid.stepY = std::ceil(grid.lenY / grid.nr_cells_y);


    return grid;
}

std::vector<std::vector<Segment_w_info>> subdivide_grid(std::vector<Segment_w_info> segments, Grid grid) {
    int curr_id = 0;
    std::vector<std::vector<Segment_w_info>> result(grid.nr_cells);

    for (int i = 0; i < segments.size(); ++i) {
        Point shortest_point;
        double shortest_dist = DBL_MAX;
        curr_id=segments[i].poly_id;
        auto it1 = segments.begin() + i;
        while (curr_id==segments[i].poly_id) {
            Point shifted_p1 = segments[i].seg.source() + Vector(grid.offset_x, grid.offset_y);
            Point shifted_p2 = segments[i].seg.target() + Vector(grid.offset_x, grid.offset_y);
            if ( CGAL::to_double ((shifted_p1 - CGAL::ORIGIN).squared_length())< shortest_dist) {
                shortest_dist = CGAL::to_double ((shifted_p1 - CGAL::ORIGIN).squared_length());
                shortest_point = shifted_p1;
            }
            else if (CGAL::to_double ((shifted_p2 - CGAL::ORIGIN).squared_length()) < shortest_dist) {
                shortest_dist = CGAL::to_double ((shifted_p2 - CGAL::ORIGIN).squared_length());
                shortest_point = shifted_p2;
            }
            i++;
            if (i>=segments.size()) break;
        }
        auto it2 = segments.begin() + i;
        int x_tile=std::floor(CGAL::to_double(shortest_point.x())/grid.stepX);
        int y_tile=std::floor(CGAL::to_double(shortest_point.y())/grid.stepY);
        int grid_nr= x_tile + y_tile * grid.nr_cells_x;
        result[grid_nr].insert(result[grid_nr].end(), it1, it2);
        i--;
    }
    return result;
}


std::vector<std::vector<Segment_w_info>> subdivide_grid_recursive(std::vector<Segment_w_info> segments, Grid grid, double power) {
    std::vector<bool> skip (grid.nr_cells, false);

    auto subdivision=subdivide_grid(segments, grid);
    for (int i = 0; i < subdivision.size(); i++) {
        if (subdivision[i].size() > pow(segments.size(), power)) {

            std::cout << "ENTERED "<<std::endl;
            auto subgrid = half_grid(subdivision[i]);
            plot_grid(subdivision[i], subgrid, "grid_"+std::to_string(i)+".svg");
            auto splitted_again = subdivide_grid(subdivision[i],subgrid);
            if (subdivision.size()>skip.size()) {
                skip.resize(subdivision.size());
            }
            skip[i]=true;
            for (auto &subset : splitted_again) {
                if (subset.empty()) continue;
                subdivision.push_back(subset);
            }
            skip.resize(subdivision.size(), false);
        }
    }

    std::vector<std::vector<Segment_w_info>> result(0);
    for (int j = 0; j < subdivision.size(); ++j) {
        if (!skip[j] && !subdivision[j].empty()) {
            //std::cout<<"Size of subdivision "<<j<<": "<<subdivision[j].size()<<std::endl;
            result.push_back(subdivision[j]);
        }
    }
    return result;
}


Grid half_grid(const std::vector<Segment_w_info> &segments) {
    Grid grid(segments);

    grid.stepX = (grid.lenX > 0.0) ? grid.lenX / 2.0 : 0.0;
    grid.stepY = (grid.lenY > 0.0) ? grid.lenY / 2.0 : 0.0;

    grid.nr_cells_x = 2;
    grid.nr_cells_y = 2;
    grid.nr_cells = 4;

    return grid;
}

void plot_grid(std::vector<Segment_w_info> segments, Grid grid, std::string filename) {
    EDGE_EXTENSION::add_outer_box(segments, 0);
    std::vector<Segment> grid_plot = EDGE_EXTENSION::filter_segments(segments);
    auto segs = grid_plot[grid_plot.size()-2];
    double p1 = CGAL::to_double(segs.source().x());
    double p2 = CGAL::to_double(segs.source().y());
    double p3 = CGAL::to_double(segs.target().x());
    double p4 = CGAL::to_double(segs.target().y());

    double minX = grid.minX;
    double minY = grid.minY;
    double maxX = grid.maxX;
    double maxY = grid.maxY;
    Point source(minX, minY+grid.stepY);
    Point target(maxX, minY+grid.stepY);
    while (source.y()<maxY) {
        grid_plot.emplace_back(source, target);
        source = source + Vector(0, grid.stepY);
        target = target + Vector(0, grid.stepY);
    }
    Point source2(minX+grid.stepX, minY);
    Point target2(minX+grid.stepX, maxY);
    while (source2.x()<maxX) {
        grid_plot.emplace_back(source2, target2);
        source2 = source2 + Vector(grid.stepX, 0);
        target2 = target2 + Vector(grid.stepX, 0);
    }
    IO_FUNCTIONS::SVG::segments_to_svg(grid_plot, filename);
}


} // namespace SUBDIVISION