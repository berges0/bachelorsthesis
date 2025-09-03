//
// Created by samuel-berge on 8/20/25.
//

#ifndef SUBDIVISION_HPP
#define SUBDIVISION_HPP
#include "utils/pch.hpp"

namespace SUBDIVISION {



std::vector<std::vector<Segment_w_info>> subdivide_random(std::vector<Segment_w_info> &segments, double subset_size);

Grid root_grid(const std::vector<Segment_w_info> &segments, double to_the_power_of);

void plot_grid(std::vector<Segment_w_info> segments, Grid grid, std::string file_name);

std::vector<std::vector<Segment_w_info>> subdivide_grid(std::vector<Segment_w_info> segments, Grid grid);

std::vector<std::vector<Segment_w_info>> subdivide_grid_recursive(std::vector<Segment_w_info> segments, Grid grid, double power);

Grid half_grid(const std::vector<Segment_w_info> &segments);

} // namespace SUBDIVISION

#endif //SUBDIVISION_HPP
