#include "io/cl_parser.hpp"

cxxopts::Options CLParser::createParser() {
    cxxopts::Options options("approx", "Approximation algorithm");

    options.add_options()
        // positional
        (INSTANCE, "Instance file name", cxxopts::value<std::string>()) //
        (OUTPUT, "Output directory", cxxopts::value<std::string>())     //
        // pseudo-optional (have default values)
        (ALPHA, "Alpha", cxxopts::value<double>()->default_value("0.1"))                  //
        (TIMELIMIT, "Maximum time in seconds", cxxopts::value<int>()->default_value("10000")) //
        (VERSION, "Algorithm version", cxxopts::value<std::string>()->default_value("0.0")) //
        (THRESHOLD, "Threshold", cxxopts::value<double>()->default_value("50")) //
        (THRESHOLDSCALE, "Threshold scale", cxxopts::value<double>()->default_value("10")) //
        (DEGREE, "Degree threshold", cxxopts::value<double>()->default_value("5")) //
        (DISTANCE, "Distance threshold", cxxopts::value<double>()->default_value("0.75"))//
        (POWER, "To the power of", cxxopts::value<double>()->default_value("0.5"))//
        (THRESHOLDVARIANT, "Threshold variant", cxxopts::value<int>()->default_value("0")) //
        // real optional (no default)
        ;
    options.parse_positional({INSTANCE, OUTPUT});
    return options;
}

CLParser::CLParser(int argc, char **argv) {
    parseResult_ = createParser().parse(argc, argv);
    if (parseResult_.count(INSTANCE) == 0 || parseResult_.count(OUTPUT) == 0) {
        throw std::invalid_argument("Missing required positional arguments");
    }
}

std::string CLParser::inputFileName() const {
    return parseResult_[INSTANCE].as<std::string>();
}

std::string CLParser::outputDirectory() const {
    return parseResult_[OUTPUT].as<std::string>();
}


std::chrono::seconds CLParser::getTimeLimit() const {
    return std::chrono::seconds(parseResult_[TIMELIMIT].as<int>());
}

double CLParser::getAlpha() const {
    return parseResult_[ALPHA].as<double>();
}

std::string CLParser::getVersion() const {
    return parseResult_[VERSION].as<std::string>();
}

double CLParser::getThreshold() const {
    return parseResult_[THRESHOLD].as<double>();
}

double CLParser::getThresholdScale() const {
    return parseResult_[THRESHOLDSCALE].as<double>();
}

double CLParser::getDegree() const {
    return parseResult_[DEGREE].as<double>();
}

double CLParser::getDistance() const {
    return parseResult_[DISTANCE].as<double>();
}

double CLParser::getPower() const {
    return parseResult_[POWER].as<double>();
}

int CLParser::getThresholdVariant() const {
    return parseResult_[THRESHOLDVARIANT].as<int>();
}