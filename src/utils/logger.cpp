//
// Created by samuel-berge on 8/11/25.
//
#include "utils/logger.hpp"


Logger::Logger(const std::string input_file, const std::string output_directory, int64_t timelimit, std::string version, double alpha) {
    alpha_=alpha;
    timer_.start();
    std::filesystem::path p(input_file);
    std::string stem = p.stem().string();
    std::string out_dir = output_directory + "/"+ stem + "_"+ version;
    std::filesystem::create_directories(out_dir);
    out_dir_stem_ = out_dir + "/" + stem + "_" + version;
    nlohmann::json j;
    std::ofstream file(out_dir_stem_+"_log.json");
    file << j.dump(4) << std::endl; // pretty print with indentation

}

std::string Logger::out_dir_stem() {
    return out_dir_stem_;
}


void Logger::start_operation() {
    start_operation_=timer_.elapsed();
}

double Logger::alpha() {
    return alpha_;
}
void Logger::end_operation(const std::string &event) {
    add(event, std::to_string(operation_duration().count()));
}

std::chrono::duration<int64_t, std::milli> Logger::operation_duration() {
    end_operation_ = timer_.elapsed();
    return std::chrono::duration_cast<std::chrono::milliseconds>(end_operation_ - start_operation_);
}

bool Logger::in_Time() {
    if (time_limit_ < timer_.elapsed().count()) {;
        add("Time limit exceeded", std::to_string(timer_.elapsed().count()));
        end();
        throw std::runtime_error("Time limit exceeded");
    }
    return (time_limit_ >= timer_.elapsed().count());
}

void Logger::end() {
    timer_.stop();
    add("Overall Time taken by algorithm: " , std::to_string(timer_.elapsed().count()));
    std::ofstream outFile(out_dir_stem_+"_log.json");
    outFile << obj_.dump(4) << std::endl;
    std::cout<< "Overall Time taken by algorithm: " << timer_.elapsed() << std::endl;
}



