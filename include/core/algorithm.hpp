//
// Created by samuel-berge on 8/20/25.
//

#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "core/arrangement.hpp"
#include "core/edge_extension.hpp"
#include "core/graph.hpp"
#include "core/max_flow.hpp"
#include "core/subdivision.hpp"
#include "io/io_functions.hpp"
#include "utils/logger.hpp"

namespace ALGORITHM{

    void run_standard(const std::string &input_filename, const std::string &output_filename, double alpha, Logger &logger);

    void run_limited(const std::string &input_filename, const std::string &output_filename, double alpha, double threshold,
        Logger &logger);

    void run_subdivision(const std::string &input_filename, const std::string &output_filename, double alpha, double subset_size,
        Logger &logger);

}



#endif //ALGORITHM_HPP
