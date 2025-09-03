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

        explicit Logger(const std::string out_file_name, int64_t timelimit);

        void add(const std::string& event, JsonPushable auto&& value) {
            obj_[event] = std::forward<decltype(value)>(value);
        }
        void start_operation();

        void end_operation(const std::string &event);

        std::chrono::duration<int64_t, std::milli> operation_duration();

        bool in_Time();

        void end();

    private:

        int64_t time_limit_ = 10000;
        std::string out_file_name_;
        Timer timer_;
        nlohmann::json obj_ = nlohmann::json::object();


        std::chrono::duration<int64_t, std::milli> start_operation_;
        std::chrono::duration<int64_t, std::milli> end_operation_;

    };


#endif //LOGGER_HPP
