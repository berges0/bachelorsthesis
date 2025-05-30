#include "edge_extension.hpp"
#include <iostream>




void EdgeExtension::buildRTree(const std::vector<SHPLoader::Point2D>& points, const std::vector<std::vector<int>>& polygons) {
    for (const auto &poly:polygons) {
        int n = poly.size();
            if (n < 2) continue;

            for (int i = 0; i < n; ++i) {
                const auto& p1 = points[poly[i]];
                const auto& p2 = points[poly[(i + 1) % n]]; // Polygon schließen

                Point pt1(p1.x, p1.y);
                Point pt2(p2.x, p2.y);
                Segment edge(pt1, pt2);
                rtree.insert(edge);
            }
    }
}




std::vector<Segment> EdgeExtension::queryNearest(double x, double y, int k) const {
    std::vector<Segment> result;

    // Create a point from the query coordinates
    Point query_point(x, y);

    // Query the R-tree for the k nearest segments
    rtree.query(
        boost::geometry::index::nearest(query_point, k),
        std::back_inserter(result)
    );

    return result;
}
