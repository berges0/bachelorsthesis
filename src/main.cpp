//
// Created by samuel-berge on 8/11/25.
//
#include <fstream>
#include <iostream>
#include <vector>
#include "io/chrono_io.hpp"
#include "io/cl_parser.hpp"
#include "io/io_functions.hpp"
#include "core/algorithm.hpp"
#include "utils/logger.hpp"

void printCommandLineFlags(int argc, char **argv, std::ostream &out = std::cout) {
    out << "Num args: " << argc << ".\nArgs: ";
    for (int i = 0; i < argc; i++) {
        out << argv[i] << ", ";
    }
    out << std::endl;
}

double quantile(const std::vector<double>& v, double p) {
    double x = p * (v.size() - 1);
    size_t i = static_cast<size_t>(x);
    double f = x - i;
    return (i + 1 < v.size())
        ? v[i] * (1 - f) + v[i + 1] * f
        : v.back();
}

// Build a bbox for a CGAL segment using your BPoint/BBox types.
static inline BBox bbox_of_segment(const Segment& s) {
    const auto& A = s.source();
    const auto& B = s.target();
    const double x1 = CGAL::to_double(A.x());
    const double y1 = CGAL::to_double(A.y());
    const double x2 = CGAL::to_double(B.x());
    const double y2 = CGAL::to_double(B.y());
    BPoint minp(std::min(x1, x2), std::min(y1, y2));
    BPoint maxp(std::max(x1, x2), std::max(y1, y2));
    return BBox(minp, maxp);
}

struct MeanStddev {
    double mean;
    double stddev;
};

// Compute mean and stddev of edge->nearest segment distances
inline MeanStddev mean_edge_to_other_polygon_distance(
    const std::vector<Segment_w_info>& segs,
    std::vector<double>* per_poly_means /*=nullptr*/,
    std::vector<double>* per_poly_stdevs /*=nullptr*/
) {
    if (segs.size() < 2) {
        throw std::runtime_error("Need at least two segments.");
    }

    // Discover number of polygons
    int max_poly_id = -1;
    for (const auto& s : segs) max_poly_id = std::max(max_poly_id, s.poly_id);
    const std::size_t poly_count = static_cast<std::size_t>(max_poly_id + 1);

    // Build R-tree
    std::vector<RItem> items; items.reserve(segs.size());
    for (int i = 0; i < static_cast<int>(segs.size()); ++i) {
        items.emplace_back(bbox_of_segment(segs[i].seg), i);
    }
    RTree rtree(items.begin(), items.end());

    // Accumulators
    double global_sum = 0.0;
    double global_sum2 = 0.0;   // sum of squares
    std::size_t global_cnt = 0;

    std::vector<double> poly_sum(poly_count, 0.0);
    std::vector<double> poly_sum2(poly_count, 0.0);
    std::vector<std::size_t> poly_cnt(poly_count, 0);

    const std::size_t maxK = 1024;

    // Main loop
    for (int i = 0; i < static_cast<int>(segs.size()); ++i) {
        const auto& srec = segs[i];
        const int my_poly = srec.poly_id;
        if (srec.seg.source() == srec.seg.target()) continue;

        double best_d2 = std::numeric_limits<double>::infinity();

        for (std::size_t k = 8; k <= maxK; k *= 2) {
            std::vector<RItem> knn;
            rtree.query(bgi::nearest(bbox_of_segment(srec.seg), k), std::back_inserter(knn));
            bool improved = false;
            for (const auto& val : knn) {
                const int j = val.second;
                if (segs[j].poly_id == my_poly) continue;
                auto d2_exact = CGAL::squared_distance(srec.seg, segs[j].seg);
                double d2 = CGAL::to_double(d2_exact);
                if (d2 < best_d2) { best_d2 = d2; improved = true; }
            }
            if (improved) break;
            if (k >= segs.size()) break;
        }
        if (!std::isfinite(best_d2)) {
            for (int j = 0; j < static_cast<int>(segs.size()); ++j) {
                if (segs[j].poly_id == my_poly) continue;
                double d2 = CGAL::to_double(CGAL::squared_distance(srec.seg, segs[j].seg));
                if (d2 < best_d2) best_d2 = d2;
            }
        }

        const double d = std::sqrt(best_d2);

        global_sum  += d;
        global_sum2 += d * d;
        ++global_cnt;

        if (my_poly >= 0 && my_poly < static_cast<int>(poly_count)) {
            poly_sum[my_poly]  += d;
            poly_sum2[my_poly] += d * d;
            ++poly_cnt[my_poly];
        }
    }

    // Means
    if (per_poly_means) {
        per_poly_means->assign(poly_count, std::numeric_limits<double>::quiet_NaN());
        for (std::size_t pid = 0; pid < poly_count; ++pid) {
            if (poly_cnt[pid]) {
                (*per_poly_means)[pid] = poly_sum[pid] / static_cast<double>(poly_cnt[pid]);
            }
        }
    }

    // Standard deviations
    if (per_poly_stdevs) {
        per_poly_stdevs->assign(poly_count, std::numeric_limits<double>::quiet_NaN());
        for (std::size_t pid = 0; pid < poly_count; ++pid) {
            if (poly_cnt[pid]) {
                double mean = poly_sum[pid] / static_cast<double>(poly_cnt[pid]);
                double var  = poly_sum2[pid] / static_cast<double>(poly_cnt[pid]) - mean * mean;
                (*per_poly_stdevs)[pid] = var > 0 ? std::sqrt(var) : 0.0;
            }
        }
    }

    if (global_cnt == 0) {
        return { std::numeric_limits<double>::quiet_NaN(),
                 std::numeric_limits<double>::quiet_NaN() };
    }

    double global_mean = global_sum / static_cast<double>(global_cnt);
    double global_var  = global_sum2 / static_cast<double>(global_cnt) - global_mean * global_mean;
    double global_std  = global_var > 0 ? std::sqrt(global_var) : 0.0;

    return { global_mean, global_std };
}


int main(int argc, char **argv) {
    printCommandLineFlags(argc, argv);
    CLParser clParser(argc, argv);

    std::cout << "Running algorithm on input file name: " << clParser.inputFileName() << '\n';
    std::cout << "Alpha is: " << clParser.getAlpha() << '\n';
    std::cout << "Timelimit is: " << clParser.getTimeLimit() << '\n';
    std::cout << "Algorithm version is: " << clParser.getVersion() << '\n';

    std::string version_input = clParser.getVersion();

    std::string version = "0";
    std::string subversion = "0";

    if (version_input.size()>1) {
        assert(version_input.size()==3);
        version = version_input.substr(0,1);
        subversion = version_input.substr(2, 1);
    }
    else if (version_input.size()==1){
        version = version_input;
    }

    Logger logger(clParser.inputFileName(), clParser.outputDirectory(),clParser.getTimeLimit().count(),
           version_input, clParser.getAlpha());

    if (version == "1") {
        Logger logger_1(clParser.inputFileName(), clParser.outputDirectory(),clParser.getTimeLimit().count(),
        version_input,clParser.getThreshold(), clParser.getThresholdScale(), clParser.getThresholdVariant());
        logger=logger_1;
    }

    if (version == "3" || version == "5") {
        Logger logger_1(clParser.inputFileName(), clParser.outputDirectory(),clParser.getTimeLimit().count(),
        version_input,clParser.getThreshold(), clParser.getThresholdScale(), clParser.getThresholdVariant(), clParser.getDegree(),
        clParser.getDistance());
        logger=logger_1;
    }

    logger.add("Input file name", clParser.inputFileName());
    logger.add("Output directory", clParser.outputDirectory());
    logger.add("Version", version_input);
    logger.add("Threshold", clParser.getThreshold());
    logger.add("Threshold scale", clParser.getThresholdScale());
    logger.add("Threshold_deg", clParser.getDegree());
    logger.add("Threshold_dist", clParser.getDistance());
    logger.add("To the power of", clParser.getPower());
    logger.add("Timelimit", clParser.getTimeLimit().count());
    logger.add("Threshold variant", clParser.getThresholdVariant());


/*    auto polygonswh = IO_FUNCTIONS::SHP::read_shp_to_pwh(clParser.inputFileName());

    CGAL::Polygon_set_2<Kernel> pset;
    pset.join(polygonswh.begin(), polygonswh.end());
    std::vector<PWH> out;
    pset.polygons_with_holes(std::back_inserter(out));

    std::cout<<"SIZE: "<< out.size()<<'\n';
    IO_FUNCTIONS::GPKG::write_to_gpkg(out,logger.out_dir_stem()+"_delaunay_combined");


    std::filesystem::path p(clParser.inputFileName());

    std::string villagename = p.stem().string();

    std::string reference_log = "/home/samuel-berge/Downloads/experiment_w_edge_elongation/output_exp1/"+ villagename + "_v0/" + villagename + "_v0_a" +
        std::to_string(clParser.getAlpha()) + "/" + villagename + "_v0_a" + std::to_string(clParser.getAlpha()) + "_log.json";

    std::ifstream in(reference_log);
    if (!in) { throw std::runtime_error("Cannot open " + reference_log); }
    nlohmann::json j;
    in >> j; // parse JSON
    in.close();

    int num_faces_arr = j["Number of faces in arrangement"];
    int num_input_polys = j["Number of polygons in input"];
    int num_new_faces=num_faces_arr-num_input_polys;
    double ref_ov = j["Objective Value"];

    nlohmann::json js;
    js["Standard_nr_faces"]=num_new_faces;
    js["OV ref_solution"]=ref_ov;

    std::vector<Segment_w_info> input_segments;

    ALGORITHM::read_in(input_segments, clParser.inputFileName(), logger);

    auto polygonswh = IO_FUNCTIONS::GPKG::read_gpkg_to_pwh(clParser.inputFileName());

    RTree rtree = SUBSTITUTE_EDGES::build_r_tree(polygonswh);

    double avg_len=EDGE_EXTENSION::compute_average_length(input_segments);

    EDGE_EXTENSION::add_outer_box(input_segments,0.01);

    std::vector<Segment_w_info> extended_segments;

    extended_segments = EDGE_EXTENSION::STANDARD::extension(input_segments);

    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(extended_segments), "after_extension.svg");

    std::cout<<"Number of segments after extension: "<<extended_segments.size()<<std::endl;
    js["Standard_nr_extended_edges"]=extended_segments.size()-4;


    auto spatially_close = PRE_PROCESS::group_by_degree_and_closeness(extended_segments, clParser.getDegree()
        ,avg_len*clParser.getDistance());

    auto spatially_close_2 = spatially_close;

    auto just_seg = EDGE_EXTENSION::filter_segments(input_segments);
    Tree tree(just_seg.begin(), just_seg.end());
    tree.build();
    tree.accelerate_distance_queries();

    SUBSTITUTE_EDGES::relink_edges(spatially_close, tree);

    auto relinked_segments = spatially_close[0];

    if (subversion=="1") {
        SUBSTITUTE_EDGES::post_prune(relinked_segments);
    }
    std::cout<<"Number of segments after relinking: "<<relinked_segments.size()<<std::endl;
    js["Relinked_nr_extended_edges"]=relinked_segments.size()-4;

    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(relinked_segments), "after_relink.svg");



    std::vector<PWH> output_data;


    logger.start_operation();
    Arrangement arr = ARRANGEMENT::build_arrangement_relinked(relinked_segments, rtree, polygonswh, logger);
    int rel_faces=arr.number_of_faces()-num_input_polys;
    js["Relinked_nr_faces_arr"]=rel_faces;
    js["Reduction of face number by relink"] =static_cast<double>(rel_faces)/num_new_faces;

    logger.start_operation();
    Graph graph = GRAPH::build_graph(arr, clParser.getAlpha(), logger);
    logger.end_operation("Building Graph (milliseconds)");
    logger.add("Number of edges in graph", std::get<0>(graph).size());


    logger.start_operation();
    std::vector<bool> max_flow_solution = MAX_FLOW::max_flow(graph, arr);
    logger.end_operation("Running Max Flow (milliseconds)");
    int nr_faces_solution = 0;
    for (int i = 0; i < max_flow_solution.size(); i++) {if (max_flow_solution[i]) nr_faces_solution++;}
    logger.add("Number of faces solution", nr_faces_solution);

    logger.start_operation();
    //IO_FUNCTIONS::all_polygons_in_solution(output_data, max_flow_solution, arr);
    auto holes_and_outer = IO_FUNCTIONS::combine_polygons(max_flow_solution, arr, logger);
    output_data = IO_FUNCTIONS::create_polygons_with_holes(holes_and_outer.first, holes_and_outer.second);
    logger.end_operation("Combining faces of solution (milliseconds)");
    logger.add("Number of polygons in solution", holes_and_outer.first.size());
    logger.add("Number of holes in solution", holes_and_outer.second.size());
    SUBSTITUTE_EDGES::connect_outer_points(spatially_close_2);



    // --- Compute area & perimeter ---
    Kernel::FT area(0);
    Kernel::FT perimeter(0);

    for (auto &pwh : output_data) {
        // outer boundary
        area += CGAL::abs(pwh.outer_boundary().area());
        for (auto eit = pwh.outer_boundary().edges_begin(); eit != pwh.outer_boundary().edges_end(); ++eit) {
            perimeter += std::sqrt(CGAL::to_double(eit->squared_length()));
        }
        // holes
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
            std::cout<<"HAS HOLE";
            area -= CGAL::abs(hit->area());
            for (auto eit = hit->edges_begin(); eit != hit->edges_end(); ++eit) {
                perimeter += std::sqrt(CGAL::to_double(eit->squared_length()));
            }
        }
    }

    double d_area = CGAL::to_double(area);
    double d_perimeter = CGAL::to_double(perimeter);

    js["Area ER"]=d_area;
    js["Perimeter ER"]=d_perimeter;
    js["Objective Value ER"]=clParser.getAlpha()*d_area + (1-clParser.getAlpha())*d_perimeter;

    double ov_det_er= clParser.getAlpha()*d_area+ (1-clParser.getAlpha())*d_perimeter;
    ov_det_er/=ref_ov;
    ov_det_er-=1;
    js["Objective det. by ER"]=ov_det_er;

    auto connected_outer_point_2 = spatially_close_2[0];


    IO_FUNCTIONS::SVG::segments_to_svg(EDGE_EXTENSION::filter_segments(connected_outer_point_2), "after_reconnect.svg");

    std::cout<<"Number of segments after reconnecting: "<<connected_outer_point_2.size()<<std::endl;


    logger.start_operation();
    Arrangement arr2 = ARRANGEMENT::build_arrangement_relinked(connected_outer_point_2, rtree, polygonswh, logger);
    int con_faces=arr2.number_of_faces()-num_input_polys;
    js["Reconnected_nr_faces_arr"]=con_faces;
    js["Reduction of face number by reconnect"] =static_cast<double>(con_faces)/num_new_faces;

    logger.start_operation();
    Graph graph2 = GRAPH::build_graph(arr2, clParser.getAlpha(), logger);
    logger.end_operation("Building Graph (milliseconds)");
    logger.add("Number of edges in graph", std::get<0>(graph2).size());


    logger.start_operation();
    std::vector<bool> max_flow_solution2 = MAX_FLOW::max_flow(graph2, arr2);
    logger.end_operation("Running Max Flow (milliseconds)");
    int nr_faces_solution2 = 0;
    for (int i = 0; i < max_flow_solution2.size(); i++) {if (max_flow_solution2[i]) nr_faces_solution2++;}
    logger.add("Number of faces solution", nr_faces_solution2);

    logger.start_operation();
    //IO_FUNCTIONS::all_polygons_in_solution(output_data, max_flow_solution, arr);
    auto holes_and_outer2 = IO_FUNCTIONS::combine_polygons(max_flow_solution2, arr2, logger);
    auto output_data2 = IO_FUNCTIONS::create_polygons_with_holes(holes_and_outer2.first, holes_and_outer2.second);
    logger.end_operation("Combining faces of solution (milliseconds)");
    logger.add("Number of polygons in solution", holes_and_outer2.first.size());
    logger.add("Number of holes in solution", holes_and_outer2.second.size());


    // --- Compute area & perimeter ---
    Kernel::FT area2(0);
    Kernel::FT perimeter2(0);

    for (auto &pwh2 : output_data2) {
        // outer boundary
        area2 += CGAL::abs(pwh2.outer_boundary().area());
        for (auto eit2 = pwh2.outer_boundary().edges_begin(); eit2 != pwh2.outer_boundary().edges_end(); ++eit2) {
            perimeter2 += std::sqrt(CGAL::to_double(eit2->squared_length()));
        }
        // holes
        for (auto hit2 = pwh2.holes_begin(); hit2 != pwh2.holes_end(); ++hit2) {
            std::cout<<"HAS HOLE";
            area2 -= CGAL::abs(hit2->area());
            for (auto eit2 = hit2->edges_begin(); eit2 != hit2->edges_end(); ++eit2) {
                perimeter2 += std::sqrt(CGAL::to_double(eit2->squared_length()));
            }
        }
    }

    double d_area2 = CGAL::to_double(area2);
    double d_perimeter2 = CGAL::to_double(perimeter2);

    js["Area CO"]=d_area2;
    js["Perimeter CO"]=d_perimeter2;
    js["Objective Value CO"]=clParser.getAlpha()*d_area2 + (1-clParser.getAlpha())*d_perimeter2;

    double ov_det_co= clParser.getAlpha()*d_area2+ (1-clParser.getAlpha())*d_perimeter2;
    ov_det_co/=ref_ov;
    ov_det_co-=1;
    js["Objective det. by CO"]=ov_det_co;


    std::filesystem::path pattt= "/home/samuel-berge/Downloads/arr_reduction/"+std::to_string(clParser.getAlpha())+"/deg" + std::to_string(clParser.getDegree()) + "_dist"+
            std::to_string(clParser.getDistance())+"/";
    if (!std::filesystem::exists(pattt)) {
        std::filesystem::create_directories(pattt); // creates all needed parent dirs
    }


    std::vector<PWH> polsys = IO_FUNCTIONS::GPKG::read_gpkg_to_pwh(clParser.inputFileName());

    std::vector<Segment_w_info> SEGS;
    IO_FUNCTIONS::pwh_to_swi(polsys,SEGS);
    //std::cout<< "Average length:" << EDGE_EXTENSION::compute_average_length(SEGS)<<std::endl;
    double avg_length=EDGE_EXTENSION::compute_average_length(SEGS);
    js["Avg_length"]=avg_length;

    std::vector<double> per_poly_means;
    std::vector<double> per_poly_stdev;
    MeanStddev mean_stdv = mean_edge_to_other_polygon_distance(SEGS, &per_poly_means, &per_poly_stdev);
    js["Mean distance"]=mean_stdv.mean;

    std::cout << "Mean: "<< mean_stdv.mean << '\n';
    std::cout << "Standard deviation from mean: "<< mean_stdv.stddev << '\n';

    double overall_area=0;

    std::vector<double> areas(0);
    for (auto poly:polsys) {
        Polygon_2 out = poly.outer_boundary();
        areas.push_back(CGAL::to_double(out.area()));
        overall_area+= CGAL::to_double(out.area());
    }
    std::sort(areas.begin(), areas.end());
    std::cout<<"MEDIAN: " << quantile(areas, 0.5);
    double avg_area=overall_area/polsys.size();
    double variance=0;
    for (auto polyy : polsys) {
        Polygon_2 outt = polyy.outer_boundary();
        variance+=std::pow((avg_area - CGAL::to_double(outt.area())),2);
    }
    variance=variance/polsys.size();




    std::cout << "Average polygon area: " << avg_area << " m\n";
    std::cout << "Standard deviation polygon area: " << std::sqrt(variance) << "m\n";

    std::ofstream out("/home/samuel-berge/Downloads/arr_reduction/"+std::to_string(clParser.getAlpha())+"/deg" + std::to_string(clParser.getDegree()) + "_dist"+
    std::to_string(clParser.getDistance())+"/"+ villagename+".json");

    out << js.dump(4); // pretty-print with 4 spaces
    out.close();

    //Per polygon (outer rings treated individually):
    for (std::size_t pid = 0; pid < per_poly_means.size(); ++pid) {
        if (std::isnan(per_poly_means[pid])) continue;
        std::cout << "poly[" << pid << "] mean: " << per_poly_means[pid] << " m\n";
    }

    std::string shp_input;
    std::string json_string;
    std::string gpkg_ouput;
    for (const auto& entry : std::filesystem::directory_iterator(clParser.inputFileName())) {
        if (entry.is_regular_file()) {
            std::filesystem::path file_path = entry.path();
            std::string ext = file_path.extension().string();

            if (ext == ".shp") {
                shp_input = file_path.parent_path().string()+"/"+file_path.stem().string();
                gpkg_ouput= file_path.parent_path().string()+"/"+file_path.stem().string()+".gpkg";
            }
            else if (ext == ".json") {
                json_string = file_path.string();
            }
        }
    }

    std::vector<PWH> polys = IO_FUNCTIONS::SHP::read_shp_to_pwh(shp_input);

    // --- Compute area & perimeter ---
    Kernel::FT area(0);
    Kernel::FT perimeter(0);

    for (auto &pwh : polys) {
        // outer boundary
        area += CGAL::abs(pwh.outer_boundary().area());
        for (auto eit = pwh.outer_boundary().edges_begin(); eit != pwh.outer_boundary().edges_end(); ++eit) {
            perimeter += std::sqrt(CGAL::to_double(eit->squared_length()));
        }
        // holes
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
            std::cout<<"HAS HOLE";
            area -= CGAL::abs(hit->area());
            for (auto eit = hit->edges_begin(); eit != hit->edges_end(); ++eit) {
                perimeter += std::sqrt(CGAL::to_double(eit->squared_length()));
            }
        }
    }

    double d_area = CGAL::to_double(area);
    double d_perimeter = CGAL::to_double(perimeter);
    std::ifstream in(json_string);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open json file " + json_string);
    }

    nlohmann::json j;
    in >> j; // parse JSON
    in.close();

    double alpha = j["Alpha"];
    j["Area"]=d_area;
    j["Perimeter"]=d_perimeter;
    j["Objective Value"]=alpha*d_area + (1-alpha)*d_perimeter;
    // Modify or add an attribute

    // Write JSON back to file (overwrite)
    std::ofstream out(json_string);
    out << j.dump(4); // pretty print with indent of 4 spaces
    out.close();

    std::cout<< " PERIMETER: " << perimeter << '\n';
    std::cout<< " AREA: " << area << '\n';


    IO_FUNCTIONS::GPKG::write_to_gpkg(polys, gpkg_ouput);
*/
    if (version == "0") {
        ALGORITHM::run_standard(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), logger);
    } else if (version == "1") {
        ALGORITHM::run_limited(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getThreshold(),
            clParser.getThresholdScale(), clParser.getThresholdVariant(), logger);
    } else if (version == "2"){
        ALGORITHM::run_preprocessed(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getDegree(),
            clParser.getDistance(), subversion, logger);
    } else if (version == "3") {
        ALGORITHM::run_edge_relink(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getThreshold(),
            clParser.getThresholdScale(), clParser.getThresholdVariant(), clParser.getDegree(), clParser.getDistance(), subversion, logger);
    } else if (version == "4") {
        ALGORITHM::run_subdivision(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(), clParser.getPower(),
            subversion, logger);
    } else if (version == "5") {
        ALGORITHM::run_outer_endpoints(clParser.inputFileName(), clParser.outputDirectory(), clParser.getAlpha(),clParser.getThreshold(),
            clParser.getThresholdScale(), clParser.getThresholdVariant(), clParser.getDegree(), clParser.getDistance(), subversion, logger);
    } else{
        std::cerr << "Unknown version: " << version << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

