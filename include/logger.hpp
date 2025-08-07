#ifndef EXECUTION_LOGGER_HPP
#define EXECUTION_LOGGER_HPP

#include <string>
#include <chrono>

class ExecutionLogger {
public:
    void startTimer();
    void stopTimer();

    void setInputPolygonCount(int count);
    void setFaceCount(int count);
    void setResultPolygonCount(int count);

    void saveToFile(const std::string& filename) const;

private:
    using Clock = std::chrono::high_resolution_clock;

    // Time
    Clock::time_point t_start, t_end;
    double durationMs = 0.0;

    // Fixed metrics
    int inputPolygons = 0;
    int triangles = 0;
    int resultPolygons = 0;

    double largestExtensionLength = 0.0;
    int kollisionen = 0;
    int flaechenErweitert = 0;
};

#endif // EXECUTION_LOGGER_HPP
