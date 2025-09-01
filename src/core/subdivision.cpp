//
// Created by samuel-berge on 8/20/25.
//

#include "core/subdivision.hpp"

#include "core/edge_extension.hpp"
#include "io/io_functions.hpp"

namespace SUBDIVISION {

std::vector<std::vector<Segment_w_info>> subdivide_random(std::vector<Segment_w_info> &segments, double subset_size) {
    int n = segments.size();
    int k = log(n)* static_cast<int>(sqrt(n));  // chunk size = sqrt(n)
    std::vector<std::vector<Segment_w_info>> result;
    for (int i=1; i<3; ++i) {
        int index = 0;
        int limit = k/i ;
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

static std::vector<double> make_coords(double minv, double maxv, double step) {
    std::vector<double> v;
    v.reserve(32);
    v.push_back(minv);
    if (step > 0) {
        const double len = maxv - minv;
        const double n_exact = len / step;        // may be non-integer
        const int n = (int)std::floor(n_exact);   // number of *full* steps that fit
        for (int i = 1; i < n; ++i) {
            v.push_back(minv + i * step);
        }
    }
    if (v.back() != maxv) v.push_back(maxv);
    return v;
}


Grid root_grid(const std::vector<Segment_w_info> &segments, double to_the_power_of) {
    Grid grid(segments);

    grid.stepX = (grid.lenX > 0.0) ? pow(grid.lenX, to_the_power_of) : 0.0;
    grid.stepY = (grid.lenY > 0.0) ? pow(grid.lenY, to_the_power_of) : 0.0;

    grid.nr_cells_x = std::ceil(grid.lenX / grid.stepX);
    grid.nr_cells_y = std::ceil(grid.lenY/ grid.stepY);
    grid.nr_cells = grid.nr_cells_x * grid.nr_cells_y;

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

void plot_grid(std::vector<Segment_w_info> segments, Grid grid) {
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
    IO_FUNCTIONS::SVG::segments_to_svg(grid_plot, "grid.svg");
}

std::vector<Segment_w_info> pwh_to_swi(std::vector<PWH> &polygons) {
    std::vector<Segment_w_info> segments;
    int seg_id=0;
    int poly_id=0;

    for (auto &poly:polygons) {
        Polygon_2 outer_boundary = poly.outer_boundary();

        int n = outer_boundary.size();
        if(n<2) continue;

        if(n==2){
            segments.emplace_back(Segment_w_info(Segment( outer_boundary[0],outer_boundary[1]), true, seg_id++, poly_id++, true, true));
            continue;
        }

        bool backward_shoot = (CGAL::orientation(outer_boundary[n - 1],outer_boundary[0],
            outer_boundary[1]) == CGAL::LEFT_TURN);

        bool forward_shoot;

        for(int i=0; i<n; i++) {

            forward_shoot = (CGAL::orientation(outer_boundary[i], outer_boundary[(i+1) % n], outer_boundary[(i+2) % n])
                == CGAL::LEFT_TURN);

            segments.emplace_back(Segment_w_info(Segment(outer_boundary[i],outer_boundary[(i+1) % n]), true, seg_id++, poly_id, backward_shoot, forward_shoot));

            backward_shoot=forward_shoot;

        }
        // Holes
        for (PWH::Hole_iterator hit = poly.holes_begin(); hit != poly.holes_end(); ++hit) {
            if (!hit->is_clockwise_oriented()) {
                hit->reverse_orientation();
                std::cout << "Had to reverse orientation of hole in polygon with ID: " << poly_id << std::endl;
            }

            const auto &hole = *hit;
            if (hole.size() < 3) continue;

            int m = hole.size();

            backward_shoot = (CGAL::orientation(hole[m - 1],hole[0],
                hole[1]) == CGAL::LEFT_TURN);

            for(int j=0; j<m; j++) {

                forward_shoot = (CGAL::orientation(hole[j], hole[(j+1) % m], hole[(j+2) % m])
                    == CGAL::LEFT_TURN);

                segments.emplace_back(Segment_w_info(Segment(hole[j],hole[(j+1) % m]), false, seg_id++,
                    poly_id, backward_shoot, forward_shoot));
                backward_shoot=forward_shoot;
            }
        }
        poly_id++;
    }
    return segments;
}

} // namespace SUBDIVISION