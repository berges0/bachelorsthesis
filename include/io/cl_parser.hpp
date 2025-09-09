//
// Created by samuel-berge on 8/11/25.
//

#ifndef CL_PARSER_HPP
#define CL_PARSER_HPP

#include <chrono>
#include <optional>

#include "cxxopts/cxxopts.hpp"

class CLParser {
public:
    CLParser(int argc, char **argv);

    std::string inputFileName() const;

    std::string outputDirectory() const;

    double getAlpha() const;

    std::chrono::seconds getTimeLimit() const;

    std::string getVersion() const;

    double getThresholdScale() const;

    double getDegree() const;

    double getDistance() const;

    double getPower() const;

    int getThresholdVariant() const;

private:
    cxxopts::ParseResult parseResult_;

    static cxxopts::Options createParser();

    inline static std::string INSTANCE = "instance";
    inline static std::string OUTPUT = "output";
    inline static std::string ALPHA = "alpha";
    inline static std::string TIMELIMIT = "timelimit";
    inline static std::string VERSION = "version";
    inline static std::string THRESHOLDSCALE = "th_scale";
    inline static std::string DEGREE = "degree";
    inline static std::string DISTANCE = "distance";
    inline static std::string POWER = "power";
    inline static std::string THRESHOLDVARIANT = "th_variant";
};

#endif //CL_PARSER_HPP
