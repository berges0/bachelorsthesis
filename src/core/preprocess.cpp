//
// Created by samuel-berge on 8/23/25.
//

#include "core/preprocess.hpp"

#include "core/edge_extension.hpp"
#include "io/io_functions.hpp"

namespace PRE_PROCESS {

std::vector<std::vector<Segment_w_info>> group_degree(const std::vector<Segment_w_info> &segments, double deg_threshold) {

    //FIRST GROUP CONTAINS SEGMENTS TO BE IGNORED

    std::vector<std::vector<Segment_w_info>> result(std::ceil(180/deg_threshold)+1);
    for (const auto &seg : segments){
        if (seg.from_poly) {
            result[0].push_back(seg);
            continue;
        }
        const auto& a = seg.seg.source();
        const auto& b = seg.seg.target();

        double dx = CGAL::to_double(b.x()) - CGAL::to_double(a.x());
        double dy = CGAL::to_double(b.y()) - CGAL::to_double(a.y());

        // atan2 gives the angle (radians) in [-π, π]
        double angle_rad = std::atan2(dy, dx);

        // convert to degrees in [0,360)
        double angle_deg = angle_rad * 180.0 / M_PI;
        if (angle_deg < 0) angle_deg += 360.0;
        // undirected...
        if (angle_deg >= 180.0) angle_deg -= 180.0;
        const int group = std::floor(angle_deg/deg_threshold);
        result[group].push_back(seg);
    }

    return result;
}



std::vector<std::vector<Segment_w_info>> spatially_close_groups(std::vector<std::vector<Segment_w_info>> &groups, double threshold_distance) {

    //FIRST GROUP CONTAINS SEGMENTS TO BE IGNORED
    std::vector<std::vector<Segment_w_info>> spatially_close_groups = {groups[0]};
    const Kernel::FT th = Kernel::FT(threshold_distance); // convert once, keep exact

    auto *skip = &groups[0];
    for (auto &group:groups) {
        if (&group==skip)continue;
        std::vector<bool> visited(group.size(), false);
        auto just_segments = EDGE_EXTENSION::filter_segments(group);
        Tree tree(just_segments.begin(), just_segments.end());
        tree.accelerate_distance_queries();

        for (int i = 0; i < group.size(); ++i) {
            if (visited[i])continue;
            visited[i] = true;

            std::vector<Segment_w_info> spatial_group(0);
            const auto &segment = group[i];
            spatial_group.push_back(segment);

            auto a = just_segments[i].source();
            auto b = just_segments[i].target();
            Kernel::FT xmin = CGAL::min(a.x(), b.x()) - th;
            Kernel::FT ymin = CGAL::min(a.y(), b.y()) - th;
            Kernel::FT xmax = CGAL::max(a.x(), b.x()) + th;
            Kernel::FT ymax = CGAL::max(a.y(), b.y()) + th;
            Kernel::Iso_rectangle_2 rect(xmin, ymin, xmax, ymax);
            std::vector<Primitive::Id> cand_ids;
            tree.all_intersected_primitives(rect, std::back_inserter(cand_ids));
            for (auto it_id : cand_ids) {
                int id = static_cast<int>(it_id - just_segments.begin());
                if (!visited[id]) {
                    const auto &neighboor = group[id];
                    if (CGAL::squared_distance(segment.seg, neighboor.seg)<th*th) {
                        visited[id] = true;
                        spatial_group.push_back(neighboor);
                    }
                }
            }

            if (spatial_group.size()>1) {
                spatially_close_groups.push_back(spatial_group);
            }
            else if (spatial_group.size()==1) {
                spatially_close_groups[0].push_back(spatial_group[0]);
            }
        }
    }
    return spatially_close_groups;
}

void shortest_wins(std::vector<std::vector<Segment_w_info>> &spatially_close_groups) {
    // sort by segment_length descending
    auto *skip = &spatially_close_groups[0];
    for (auto &vec : spatially_close_groups) {
        if (&vec==skip)continue;
        std::sort(vec.begin(), vec.end(),
            [](const Segment_w_info &a, const Segment_w_info &b) {
                return a.seg.squared_length() < b.seg.squared_length();
            }
        );
    }
    for (int g = 1; g<spatially_close_groups.size(); ++g) {
        assert (spatially_close_groups[g].size()>0);
        spatially_close_groups[0].push_back(spatially_close_groups[g][0]);
    }
}

void longest_wins(std::vector<std::vector<Segment_w_info>> &spatially_close_groups) {
    // sort by segment_length descending
    auto *skip = &spatially_close_groups[0];
    for (auto &vec : spatially_close_groups) {
        if (&vec==skip)continue;
        std::sort(vec.begin(), vec.end(),
            [](const Segment_w_info &a, const Segment_w_info &b) {
                return a.seg.squared_length() > b.seg.squared_length();
            }
        );
    }
    for (int g = 1; g<spatially_close_groups.size(); ++g) {
        assert (spatially_close_groups[g].size()>0);
        spatially_close_groups[0].push_back(spatially_close_groups[g][0]);
    }
}

void longest_and_shortest_wins(std::vector<std::vector<Segment_w_info>> &spatially_close_groups) {
    // sort by segment_length descending
    auto *skip = &spatially_close_groups[0];
    for (auto &vec : spatially_close_groups) {
        if (&vec==skip)continue;
        std::sort(vec.begin(), vec.end(),
            [](const Segment_w_info &a, const Segment_w_info &b) {
                return a.seg.squared_length() > b.seg.squared_length();
            }
        );
    }
    for (int g = 1; g<spatially_close_groups.size(); ++g) {
        spatially_close_groups[0].push_back(spatially_close_groups[g][0]);
        if (spatially_close_groups[g].size()>1) {
            spatially_close_groups[0].push_back(spatially_close_groups[g].back());
        }
    }
}

void longest_mid_shortest_wins(std::vector<std::vector<Segment_w_info>> &spatially_close_groups) {
    // sort by segment_length descending
    auto *skip = &spatially_close_groups[0];
    for (auto &vec : spatially_close_groups) {
        if (&vec==skip)continue;
        std::sort(vec.begin(), vec.end(),
            [](const Segment_w_info &a, const Segment_w_info &b) {
                return a.seg.squared_length() > b.seg.squared_length();
            }
        );
    }
    for (int g = 1; g<spatially_close_groups.size(); ++g) {
        spatially_close_groups[0].push_back(spatially_close_groups[g][0]);
        if (spatially_close_groups[g].size()>2) {
            spatially_close_groups[0].push_back(spatially_close_groups[g][spatially_close_groups[g].size()/2]);
        }
        if (spatially_close_groups[g].size()>1) {
            spatially_close_groups[0].push_back(spatially_close_groups[g].back());
        }
    }
}
}





