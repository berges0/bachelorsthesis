#include "algorithm.hpp"


void algorithm::run(std::string input_path) {

    // Load instance using the library
    auto data = SHPLoader::ReadShapeFileToPoint2D(input_path);
    std::cout << "Successfully loaded instance from: " << input_path << std::endl;
    std::cout << "Points: " << data.first.size() << std::endl;
    std::cout << "Polygons: " << data.second.size() << std::endl;
    auto transformed_to_segments = build_segment_vector(data.first, data.second);
    std::vector<Segment_w_info> segments_w_info = transformed_to_segments.first;
    add_box(segments_w_info, transformed_to_segments.second);
    std::vector<Segment> segments= segs_wo_info(segments_w_info);
    /*
    std::vector<Segment> segments;
    segments.emplace_back(Point(1.0,0.0), Point(1.0, 4.0));
    segments.emplace_back(Point(4.0, 0.0), Point(4.0, 4.0));
    segments.emplace_back(Point(2.0, 1.0), Point(3.0, 1.0));
    segments.emplace_back(Point(2.0, 3.0), Point(3.0, 3.0));
    segments.emplace_back(Point(2.0, 0.0), Point(3.0, 0.0));
    segments.emplace_back(Point(2.0, 4.0), Point(3.0, 4.0));
    segments.emplace_back(Point(3.0, 1.0), Point(3.0, 3.0));
    std::vector<Segment_w_info> segments_w_info;
    int id=0;
    for (auto seg : segments) {
        segments_w_info.emplace_back(Segment_w_info(seg, true, id));
        ++id;
    }*/

    segments_to_svg(segments,"original.svg");
    auto extended_segs = ray_shoot_intersection(segments_w_info);
    std::cout<<"Extended Segments size: "<<extended_segs.size()<<std::endl;
    segments_to_svg(segs_wo_info(extended_segs),"extended.svg");

    Arrangement arr = build_arrangement(extended_segs);
    auto graph = build_graph(arr, 0.001);
    auto output_data = max_flow(graph,arr);
    segments_to_svg(output_data, "output_segs.svg");
    //write_output_polygons(output_data, out_path);
    std::cout << "Arrangement has "
          << arr.number_of_vertices() << " vertices, "
          << arr.number_of_edges() << " edges, "
          << arr.number_of_faces() << " faces." << std::endl;
}