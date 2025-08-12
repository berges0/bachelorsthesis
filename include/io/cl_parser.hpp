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

    std::string outputFileName() const;

    double getAlpha() const;

    std::chrono::seconds getTimeLimit() const;

    std::string getVersion() const;

private:
    cxxopts::ParseResult parseResult_;

    static cxxopts::Options createParser();

    inline static std::string INSTANCE = "instance";
    inline static std::string OUTPUT = "output";
    inline static std::string ALPHA = "alpha";
    inline static std::string TIMELIMIT = "timelimit";
    inline static std::string VERSION = "version";
};

#endif //CL_PARSER_HPP
