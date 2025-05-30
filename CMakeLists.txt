cmake_minimum_required(VERSION 3.20)
project(edge_aligned_aggregation)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add external lib
add_subdirectory(lib/shapefileaggregationloader)
# Add Boost package
find_package(Boost REQUIRED)


# Add your own code
add_executable(aggregator src/main.cpp
        src/edge_extension.cpp
        src/edge_extension.cpp)

target_include_directories(aggregator PRIVATE include)

target_include_directories(aggregator PRIVATE ${Boost_INCLUDE_DIRS})

# Link your code with the library
target_link_libraries(aggregator
        PRIVATE instance_loader_lib
)
