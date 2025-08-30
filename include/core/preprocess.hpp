//
// Created by samuel-berge on 8/23/25.
//

#ifndef PREPROCESS_HPP
#define PREPROCESS_HPP
#include <vector>

#include "utils/pch.hpp"

namespace PRE_PROCESS {


std::vector<std::vector<Segment_w_info>> group_degree(const std::vector<Segment_w_info> &segments, double deg_threshold);

void shortest_wins(std::vector<std::vector<Segment_w_info>> &spatially_close_groups);

void longest_wins(std::vector<std::vector<Segment_w_info>> &spatially_close_groups);

void longest_mid_shortest_wins(std::vector<std::vector<Segment_w_info>> &spatially_close_groups);

void longest_and_shortest_wins(std::vector<std::vector<Segment_w_info>> &spatially_close_groups);

std::vector<std::vector<Segment_w_info>> spatially_close_groups(std::vector<std::vector<Segment_w_info>> &groups,
    double threshold_distance);

std::vector<Segment_w_info> merge_spatially_close(std::vector<Segment_w_info> &input_segments,
    std::vector<std::vector<Segment_w_info>> &groups, double threshold_distance);

std::vector<Segment_w_info> merge_edges(std::vector<std::vector<Segment_w_info>> spatial_groups,
    std::vector<Segment_w_info> &input_segments, std::vector<Segment_w_info> &unmodified);
}



#endif //PREPROCESS_HPP
