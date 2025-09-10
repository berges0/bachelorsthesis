//
// Created by samuel-berge on 8/10/25.
//

#ifndef LOGGER_HPP
#define LOGGER_HPP
#include "io/chrono_io.hpp"
#include "utils/json.hpp"
#include <fstream>

template<class T>
concept JsonPushable =
requires (const T& v) { json(v); } ||
requires (nlohmann::json& j, const T& v) { to_json(j, v); };

    class Logger {
    public:

        explicit Logger(const std::string input_file, const std::string output_directory, int64_t timelimit, std::string version);

        Logger(const std::string input_file, const std::string output_directory, int64_t timelimit, std::string version,
            double threshold, double th_scale, int threshold_variant);

        Logger(Logger &logger, double alpha);


        void add(const std::string& event, JsonPushable auto&& value) {
            obj_[event] = std::forward<decltype(value)>(value);
        }

        std::string out_dir_stem();

        void start_operation();

        void end_operation(const std::string &event);

        std::chrono::duration<int64_t, std::milli> operation_duration();

        bool in_Time();

        void end();

        nlohmann::json obj();


        int64_t time();

        void stop_time();

        void start();

        void end_sub_log();

        double alpha();


    private:
        std::string input_file;
        std::string out_dir_stem_;

        int64_t time_limit_ = 10000;
        int64_t until_now;
        Timer timer_;
        nlohmann::json obj_ = nlohmann::json::object();
        double alpha_;

        std::chrono::duration<int64_t, std::milli> start_operation_;
        std::chrono::duration<int64_t, std::milli> end_operation_;

    };


#endif //LOGGER_HPP
