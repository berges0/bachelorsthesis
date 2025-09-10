//
// Created by samuel-berge on 8/11/25.
//
#include "utils/logger.hpp"


Logger::Logger(const std::string input_file, const std::string output_directory, int64_t timelimit, std::string version) {
    timer_.start();
    std::filesystem::path p(input_file);
    std::string stem = p.stem().string();
    std::string out_dir = output_directory + "/"+ stem + "_v"+ version;
    std::filesystem::create_directories(out_dir);
    out_dir_stem_ = out_dir + "/" + stem + "_v" + version;
}

Logger::Logger(const std::string input_file, const std::string output_directory, int64_t timelimit, std::string version,
    double threshold, double th_scale, int threshold_variant) {
    std::string suffix;
    if (threshold_variant ==0) {
        suffix = "_tv0_t" + std::to_string(threshold);
    }
    else if (threshold_variant == 1) {
        suffix = "_tv1_ts" + std::to_string(th_scale);
    }
    else if (threshold_variant == 2) {
        suffix = "_tv2_t" + std::to_string(threshold) + "_ts" + std::to_string(th_scale);
    }
    else{
        throw std::runtime_error("threshold variant not recognized");
    }
    timer_.start();
    std::filesystem::path p(input_file);
    std::string stem = p.stem().string();
    std::string out_dir = output_directory + "/"+ stem + "_v"+ version+ suffix;
    std::filesystem::create_directories(out_dir);
    out_dir_stem_ = out_dir + "/" + stem + "_v" + version+ suffix;
}

std::string Logger::out_dir_stem() {
    return out_dir_stem_;
}


void Logger::start_operation() {
    start_operation_=timer_.elapsed();
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

nlohmann::json Logger::obj() {
    return obj_;
}

void Logger::stop_time() {
    timer_.stop();
}

int64_t Logger::time() {
    return timer_.elapsed().count();
}

Logger::Logger(Logger &logger, double alpha) {
    alpha_ = alpha;
    std::string dir = logger.out_dir_stem()+"_a" +std::to_string(alpha);
    std::filesystem::create_directories(dir);
    std::filesystem::path p(dir);
    std::string stem = p.filename().string();
    out_dir_stem_ = dir + "/" + stem;
    obj_ = logger.obj();
    until_now = logger.time();
}


void Logger::start() {
    timer_.start();
}

double Logger::alpha() {
    return alpha_;
}

void Logger::end_sub_log() {
    timer_.stop();
    int64_t overall=until_now + timer_.elapsed().count();
    add("Overall Time taken by algorithm: " , std::to_string(overall));
    std::ofstream outFile(out_dir_stem_+"_log.json");
    outFile << obj_.dump(4) << std::endl;
    std::cout<< "Overall Time taken by algorithm: " << std::to_string(overall) << std::endl;
}




