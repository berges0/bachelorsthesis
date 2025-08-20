//
// Created by samuel-berge on 8/10/25.
//
#ifndef IO_FUNCTIONS_HPP
#define IO_FUNCTIONS_HPP

#include "utils/pch.hpp"

namespace IO_FUNCTIONS {

//CONTROVERSION FROM SHAPEFILEAGGREGATIONLOADER TO CGAL DATASTRUCTURE

    std::vector<Segment_w_info> read_in_segment(const std::vector<SHPLoader::Point2D>& points,
        const std::vector<std::vector<int>>& polygons);


//SVG FILES

    void segments_to_svg(const std::vector<Segment>& segments, const std::string& filename);

    void polygons_to_svg(const std::vector<Polygon_with_holes_2>& polygons, const std::string& filename);


//SHP FILES

    void write_to_shp(const std::vector<Polygon_with_holes_2>& polygons, const std::string& filename);

    void writeToShapeFile(std::vector<Polygon_with_holes_2> polys, std::string path);

//CONVERT SOLUTION POLYGONS TO CONTIGUOUS POLYGONS

    const std::pair <std::vector<Polygon_2>, std::vector<Polygon_2>>  combine_polygons(const std::vector<bool>
        &in_solution, Arrangement &arr);

    const Polygon_2 get_contiguous_boundary(Arrangement::Halfedge_handle &edge, const std::vector<bool> &in_solution);

    const std::map<int, std::vector<Polygon_2>> locate_holes(const std::vector<Polygon_2> &outer_boundaries,
        const std::vector<Polygon_2> &holes);

    const std::vector<Polygon_with_holes_2> create_polygons_with_holes(const std::vector<Polygon_2> &outer_boundaries,
        const std::vector<Polygon_2> &holes);

    const std::vector<Polygon_with_holes_2> cgal_combines(const std::vector<Polygon_2> &polygons);

}


#endif //IO_FUNCTIONS_HPP
