#include "algorithm.hpp"

int main(int argc, char* argv[]) {
    std::cout << "Hello from Edge-Aligned Aggregation!" << std::endl;

    // Require at least 1 argument: input file path
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_shapefile_path> [alpha]" << std::endl;
        return 1;
    }

    std::string input_path = argv[1];

    // Check if input file exists
    if (!std::filesystem::exists(input_path)) {
        std::cerr << "Error: Input file does not exist: " << input_path << std::endl;
        return 2;
    }

    // Default alpha value
    double alpha = 0.0003;

    // Parse optional alpha argument
    if (argc >= 3) {
        try {
            alpha = std::stod(argv[2]);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid alpha value. Must be a number. Received: " << argv[2] << std::endl;
            return 3;
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: Alpha value out of range. Received: " << argv[2] << std::endl;
            return 4;
        }
    }

    std::cout << "Alpha is: " << alpha << std::endl;

    algorithm algo;
    algo.run(input_path, alpha);
    std::string output_path = "/home/samuel-berge/Work/bachelorsthesis/build/shp/shp_output.shp";
    std::string command1 = "qgis " + input_path + " &";
    std::string command2 = "qgis " + output_path + " &";
    int i = std::system(command1.c_str());
    int j = std::system(command2.c_str());

    return 0;
}
