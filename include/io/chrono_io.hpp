//
// Created by samuel-berge on 8/11/25.
//

#ifndef CHRONO_IO_HPP
#define CHRONO_IO_HPP

#include <chrono>
#include <concepts>
#include <iostream>
#include <type_traits>

template <typename T>
    requires std::same_as<T, std::nano>
inline std::string suffix() {
    return "ns";
}

template <typename T>
    requires std::same_as<T, std::micro>
inline std::string suffix() {
    return "µs";
}

template <typename Ratio>
    requires std::same_as<Ratio, std::milli>
inline std::string suffix() {
    return "ms";
}

// no ratio means seconds
template <typename Rep>
std::stringstream &operator<<(std::stringstream &os, const std::chrono::duration<Rep> &duration) {
    os << duration.count() << "s";
    return os;
}

template <typename Rep>
std::ostream &operator<<(std::ostream &os, const std::chrono::duration<Rep> &duration) {
    os << duration.count() << "s";
    return os;
}

// ratio allows to deduce unit
template <typename Rep, typename Ratio>
std::stringstream &operator<<(std::stringstream &os, const std::chrono::duration<Rep, Ratio> &duration) {
    os << duration.count() << suffix<Ratio>();
    return os;
}

template <typename Rep, typename Ratio>
std::ostream &operator<<(std::ostream &os, const std::chrono::duration<Rep, Ratio> &duration) {
    os << duration.count() << suffix<Ratio>();
    return os;
}

/**
 * A timer for running time measurements.
 * Based on the one NetworKit ships with.
 */
class Timer {
public:
    Timer() = default;

    /**
     * Start the clock.
     * Returns the time at which the instance was started.
     */
    std::chrono::steady_clock::time_point start() noexcept;

    /**
     * Stops the clock permanently for the instance of the Timer.
     * Returns the time at which the instance was stopped.
     */
    std::chrono::steady_clock::time_point stop() noexcept;

    /**
     * Returns a chrono::duration since the Timer was started.
     * If stop() was called, the duration is between the start() and stop()
     * calls is returned.
     */
    std::chrono::duration<int64_t, std::milli> elapsed() const noexcept;

    /**
     * Returns the time at which the instance was started.
     */
    std::chrono::steady_clock::time_point startTime() const noexcept;

    /**
     * Returns the time at which the instance was stopped.
     */
    std::chrono::steady_clock::time_point stopTime() const noexcept;

protected:
    bool running_{false};                           //!< true if timer has been started and not stopped after that
    std::chrono::steady_clock::time_point started_; //!< time at which timer has been started
    std::chrono::steady_clock::time_point stopped_; //!< time at which timer has been stopped

    /// If running returns now, otherwise the stop time
    std::chrono::steady_clock::time_point stopTimeOrNow() const noexcept;
};

#endif //CHRONO_IO_HPP
