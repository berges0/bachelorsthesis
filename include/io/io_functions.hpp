//
// Created by samuel-berge on 8/10/25.
//
#ifndef IO_FUNCTIONS_HPP
#define IO_FUNCTIONS_HPP

#include "gdal/gdal_priv.h"
#include "gdal/ogrsf_frmts.h"
#include "utils/logger.hpp"
#include "utils/pch.hpp"
#include "utils/logger.hpp"

namespace IO_FUNCTIONS {

    // SHP I/O
    namespace SHP {

        void read(const std::string filename, std::vector<Segment_w_info>& segments, Logger &logger);
        void write_to_shp(std::vector<PWH> polys, std::string path);

    } // namespace SHP

    // GPKG I/O
    namespace GPKG {

        void read(std::string filename, std::vector<Segment_w_info>& input_segments, Logger &logger);

        std::vector<PWH> read_gpkg_to_pwh(std::string filename);

        void write_to_gpkg(const std::vector<PWH>& polys, const std::string& path);

        // Helper: OGR -> CGAL
        PWH OGRPolygonToCGAL(OGRPolygon* poPolygon);

    } // namespace GPKG

    // CGAL PWH -> Segment_w_info
    void pwh_to_swi(std::vector<PWH>& polygons, std::vector<Segment_w_info>& segments);

    // SVG
    namespace SVG {
        void segments_to_svg(const std::vector<Segment>& segments, const std::string& filename);

        void polygons_to_svg(const std::vector<PWH>& polygons, const std::string& filename);

    }// namespace SVG

    void all_polygons_in_solution(std::vector<PWH> &    all_polys_in_sol, std::vector<bool> max_flow_solution, const Arrangement& arr);

    void arrangement_as_polys(std::vector<PWH> &all_polys, const Arrangement& arr);

    // Polygon-Komposition
    const std::pair<std::vector<Polygon_2>, std::vector<Polygon_2>> combine_polygons(const std::vector<bool>& in_solution,
        Arrangement& arr);

    const Polygon_2 get_contiguous_boundary(Arrangement::Halfedge_handle& edge, const std::vector<bool>& in_solution);

    const std::map<int, std::vector<Polygon_2>> locate_holes(const std::vector<Polygon_2>& outer_boundaries,
        const std::vector<Polygon_2>& holes);

    const std::vector<PWH> create_polygons_with_holes(const std::vector<Polygon_2>& outer_boundaries,
                               const std::vector<Polygon_2>& holes);

    const std::vector<PWH> cgal_combines(const std::vector<Polygon_2>& polygons);

} // namespace IO_FUNCTIONS


#endif //IO_FUNCTIONS_HPP
