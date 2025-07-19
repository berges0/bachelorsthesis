#include "algorithm.hpp"


void algorithm::run(std::string input_path, double alpha) {

    // Load instance using the library
    auto data = SHPLoader::ReadShapeFileToPoint2D(input_path);
    std::cout << "Successfully loaded instance from: " << input_path << std::endl;
    std::cout << "Points: " << data.first.size() << std::endl;
    std::cout << "Polygons: " << data.second.size() << std::endl;
    auto transformed_to_segments = build_segment_vector(data.first, data.second);

    //auto transformed_to_segments = build_segment_vector(points, polygons);
    std::vector<Segment_w_info> segments_w_info = transformed_to_segments.first;
    add_box(segments_w_info, transformed_to_segments.second);
    std::vector<Segment> segments= segs_wo_info(segments_w_info);
    segments_to_svg(segments,"original.svg");
    auto extended_segs = ray_shoot_intersection(segments_w_info);
    std::cout<<"Extended Segments size: "<<extended_segs.size()<<std::endl;
    segments_to_svg(segs_wo_info(extended_segs),"extended.svg");

    Arrangement arr = build_arrangement(extended_segs);
    auto graph = build_graph(arr, alpha);
    auto max_sol = max_flow(graph,arr);
    auto output_data = combined_output_polygons(max_sol, arr);
    //auto output_data = output_polygons(max_sol, arr);
    polygons_to_svg(output_data, "output_polys.svg");
    write_to_shp(output_data, "/home/berges0/Uni/bachelorsthesis/build/shp/shp_output.shp");
    //write_output_polygons(output_data, out_path);
    std::cout << "Arrangement has "
          << arr.number_of_vertices() << " vertices, "
          << arr.number_of_edges() << " edges, "
          << arr.number_of_faces() << " faces." << std::endl;
}