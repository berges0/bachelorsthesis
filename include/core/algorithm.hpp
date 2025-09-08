//
// Created by samuel-berge on 8/20/25.
//

#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "core/arrangement.hpp"
#include "core/edge_extension.hpp"
#include "core/graph.hpp"
#include "core/max_flow.hpp"
#include "core/preprocess.hpp"
#include "core/subdivision.hpp"
#include "core/substitute_edges.hpp"
#include "io/io_functions.hpp"
#include "utils/logger.hpp"

namespace ALGORITHM{

    void read_in (std::vector<Segment_w_info> &input_segments, const std::string &input_filename, Logger &logger);

    void aggregate(std::vector<PWH> &output_data, const std::vector<Segment_w_info> &input, double alpha, Logger &logger);

    void run_standard(const std::string &input_filename, const std::string &output_filename, double alpha, Logger &logger);

    void run_limited(const std::string &input_filename, const std::string &output_filename, double alpha, double threshold_scale,
        int th_variant, Logger &logger);

    void run_subdivision(const std::string &input_filename, const std::string &output_filename, double alpha, double to_the_power_of,
        std::string subversion, Logger &logger);

    void run_preprocessed(const std::string &input_filename, const std::string &output_filename, double alpha, double degree,
        double distance, std::string subversion, Logger &logger);

    void run_edge_relink(const std::string &input_filename, const std::string &output_filename, double alpha, double threshold_scale,
        int th_variant, double degree, double distance, std::string subversion, Logger &logger);

    void run_outer_endpoints(const std::string &input_filename, const std::string &output_filename, double alpha, double threshold_scale,
        int th_variant, double degree, double distance, std::string subversion, Logger &logger);

}



#endif //ALGORITHM_HPP
