cmake_minimum_required(VERSION 3.20)
project(edge_aligned_aggregation)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add external subdirectory (your shapefile loader lib)
add_subdirectory(lib/shapefileaggregationloader)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)

# Find Boost
find_package(Boost REQUIRED)

# Find CGAL
find_package(CGAL REQUIRED)

# Find GDAL
find_package(GDAL REQUIRED)

# Your executable
add_executable(aggregator
        src/main.cpp
        src/core/arrangement.cpp
        src/core/edge_extension.cpp
        src/core/graph.cpp
        src/core/max_flow.cpp
        src/io/chrono_io.cpp
        src/io/cl_parser.cpp
        src/io/io_functions.cpp
        src/utils/logger.cpp
        src/core/algorithm.cpp
        include/core/subdivision.hpp
        src/core/subdivision.cpp
        include/core/preprocess.hpp
        src/core/preprocess.cpp
        include/core/substitute_edges.hpp
        src/core/substitute_edges.cpp
        src/utils/metric_simplification.cpp
)
# Precompile headers
target_precompile_headers(aggregator PRIVATE include/utils/pch.hpp)

# Include paths
target_include_directories(aggregator PRIVATE
        include
        ${GDAL_INCLUDE_DIRS}   # GDAL include directories
        ${Boost_INCLUDE_DIRS}
        ${CGAL_INCLUDE_DIRS}  # CGAL include directories (should be found automatically)
)

# Link libraries
target_link_libraries(aggregator
        PRIVATE GDAL::GDAL  # Link GDAL libraries
        PRIVATE CGAL::CGAL
        PRIVATE instance_loader_lib  # Ensure this is linked properly from the submodule
)
