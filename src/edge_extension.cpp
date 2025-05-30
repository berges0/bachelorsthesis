#include "edge_extension.hpp"
#include <iostream>




void EdgeExtension::buildRTree(const std::vector<SHPLoader::Point2D>& points, const std::vector<std::vector<int>>& polygons) {
    for (size_t p = 0; p < polygons.size(); p++) {
        int n = polygons[p].size();
            if (n < 2) continue;

            for (int i = 0; i < n; ++i) {
                const auto& p1 = points[polygons[p][i]];
                const auto& p2 = points[polygons[p][(i + 1) % n]]; // Polygon schließen

                Point pt1(p1.x, p1.y);
                Point pt2(p2.x, p2.y);
                MySegment edge(pt1, pt2, {p});
                rtree.insert(edge);
            }
    }
}


std::vector<MySegment> EdgeExtension::queryNearest(double x, double y, int k) const {
    std::vector<MySegment> result;

    // Create a point from the query coordinates
    Point query_point(x, y);

    // Query the R-tree for the k nearest segments
    rtree.query(
        boost::geometry::index::nearest(query_point, k),
        std::back_inserter(result)
    );
    int count=0;
    for (const auto& segment : rtree) {
        std::cout << "Points (" << segment.geometry.first.x() <<", " <<  segment.geometry.first.y() << ") and (" << segment.geometry.second.x()<< ", " <<  segment.geometry.second.y()<< ") from Polygon " << segment.polygon_ids[0]<< std::endl;
        if (count == 100) {
            break;
        }
        count++;

    }


    return result;
}
