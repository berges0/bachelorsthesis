//
// Created by samuel-berge on 8/11/25.
//
#include "utils/logger.hpp"

Logger::Logger(const std::string out_file_name, int64_t timelimit):out_file_name_(out_file_name), time_limit_(timelimit) {
    timer_.start();
    nlohmann::json j;
    std::ofstream file(out_file_name_);
    file << j.dump(4) << std::endl; // pretty print with indentation
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
    std::ofstream outFile(out_file_name_);
    outFile << obj_.dump(4) << std::endl;
    std::cout<< "Overall Time taken by algorithm: " << timer_.elapsed() << std::endl;
}



