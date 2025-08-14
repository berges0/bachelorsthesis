#include "io/chrono_io.hpp"

std::chrono::steady_clock::time_point Timer::start() noexcept {
    started_ = std::chrono::steady_clock::now();
    running_ = true;
    return started_;
}

std::chrono::steady_clock::time_point Timer::stop() noexcept {
    stopped_ = std::chrono::steady_clock::now();
    running_ = false;
    return stopped_;
}

std::chrono::duration<int64_t, std::milli> Timer::elapsed() const noexcept {
    auto diff = stopTimeOrNow() - started_;
    return std::chrono::duration_cast<std::chrono::duration<int64_t, std::milli>>(diff);
}

std::chrono::steady_clock::time_point Timer::startTime() const noexcept {
    return started_;
}

std::chrono::steady_clock::time_point Timer::stopTime() const noexcept {
    return stopped_;
}

std::chrono::steady_clock::time_point Timer::stopTimeOrNow() const noexcept {
    return running_ ? std::chrono::steady_clock::now() : stopped_;
}