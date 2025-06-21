#include "algorithm.hpp"


int main(int argc, char* argv[]) {
    std::cout << "Hello from Edge-Aligned Aggregation!" << std::endl;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_shapefile_path>" << std::endl;
        return 1;
    }

    std::string input_path = argv[1];

    // Check if file exists
    if (!std::filesystem::exists(input_path)) {
        std::cerr << "Error: Input file does not exist: " << input_path << std::endl;
        return 2;
    }

    algorithm algo;
    algo.run(input_path);

    return 0;
}