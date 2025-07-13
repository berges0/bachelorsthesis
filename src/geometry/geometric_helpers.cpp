
#include "algorithm.hpp"

double algorithm::calculate_length(const Arrangement::Halfedge_const_handle& he) {
    auto p1 = he->source()->point();
    auto p2 = he->target()->point();

    auto dx = p2.x() - p1.x();
    auto dy = p2.y() - p1.y();

    auto squared_length = dx*dx + dy*dy;
    return std::sqrt(CGAL::to_double(squared_length));
}