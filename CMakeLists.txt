cmake_minimum_required(VERSION 3.20)
project(edge_aligned_aggregation)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add external subdirectory (your shapefile loader lib)
add_subdirectory(lib/shapefileaggregationloader)

# Find Boost
find_package(Boost REQUIRED)

# Find CGAL
find_package(CGAL REQUIRED)

# Your executable
add_executable(aggregator
        src/main.cpp
        src/algorithm.cpp
        src/algo_geometry.cpp
        src/arrangement.cpp
        src/graphbuild.cpp
        src/io_functions.cpp
)
target_precompile_headers(aggregator PRIVATE include/pch.hpp)
# Include paths
target_include_directories(aggregator PRIVATE
        include
        ${Boost_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(aggregator
        PRIVATE CGAL::CGAL
        PRIVATE instance_loader_lib
)
