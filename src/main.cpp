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

    if (version == "1" || subversion =="1") {
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
    std::vector<PWH> polsys = IO_FUNCTIONS::GPKG::read_gpkg_to_pwh(clParser.inputFileName());

    std::vector<Segment_w_info> SEGS;
    IO_FUNCTIONS::pwh_to_swi(polsys,SEGS);
    std::cout<< "Average length:" << EDGE_EXTENSION::compute_average_length(SEGS)<<std::endl;
    /*

    std::vector<double> per_poly_means;
    std::vector<double> per_poly_stdev;
    MeanStddev mean_stdv = mean_edge_to_other_polygon_distance(SEGS, &per_poly_means, &per_poly_stdev);
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

