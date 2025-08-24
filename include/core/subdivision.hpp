//
// Created by samuel-berge on 8/20/25.
//

#ifndef SUBDIVISION_HPP
#define SUBDIVISION_HPP
#include "utils/pch.hpp"

namespace SUBDIVISION {



std::vector<std::vector<Segment_w_info>> subdivide_random(std::vector<Segment_w_info> &segments, double subset_size);



Grid root_grid(const std::vector<Segment_w_info> &segments, double to_the_power_of);

Grid log_grid(const std::vector<Segment_w_info> &segments, double log_base);

Grid root_log_grid(const std::vector<Segment_w_info> &segments, double to_the_power_of, double log_base);

void plot_grid(std::vector<Segment_w_info> segments, Grid grid);

std::vector<std::vector<Segment_w_info>> subdivide_grid(std::vector<Segment_w_info> segments, int nr_polys, Grid grid);

std::vector<Segment_w_info> pwh_to_swi(std::vector<Polygon_with_holes_2> &polygons);

} // namespace SUBDIVISION

#endif //SUBDIVISION_HPP
