//
// Created by samuel-berge on 8/23/25.
//

#include "core/preprocess.hpp"

#include "core/edge_extension.hpp"
#include "io/io_functions.hpp"

namespace PRE_PROCESS {
/*
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
*/
std::vector<std::vector<Segment_w_info>> group_by_degree_and_closeness(std::vector<Segment_w_info> &segments, double degree,
    double distance) {
    std::vector<Segment_w_info> input_segments;
    std::vector<Segment_w_info> extended_segments;
    std::vector<std::vector<Segment_w_info>> result(1);

    for (const auto &seg : segments) {
        if (seg.from_poly || seg.poly_id==-4) { //-4 is for outer box
            result[0].push_back(seg);
        }
        else {
            extended_segments.push_back(seg);
        }
    }

    double sqd_distance = distance * distance;
    double theta_max = degree * M_PI / 180.0; // convert to radians
    // Build boxes + R-tree
    std::vector<RItem> items; items.reserve(extended_segments.size());
    std::vector<BBox>  boxes(extended_segments.size());
    for (int i = 0; i < (int)extended_segments.size(); ++i) {
        auto bb = extended_segments[i].seg.bbox();
        boxes[i] = BBox(BPoint(bb.xmin(), bb.ymin()), BPoint(bb.xmax(), bb.ymax()));
        items.emplace_back(boxes[i], i);
    }
    RTree rtree(items.begin(), items.end()); // bulk load

    // Precompute angles in [0, pi) so 180° ≡ 0°
    std::vector<double> ang(extended_segments.size());
    for (int i = 0; i < extended_segments.size(); ++i) {
        auto v = extended_segments[i].seg.to_vector();
        double a = std::atan2(CGAL::to_double(v.y()), CGAL::to_double(v.x()));
        if (a < 0) a += M_PI;        // normalize to [0, pi)
        if (a >= M_PI) a -= M_PI;
        ang[i] = a;
    }

    // Undirected adjacency
    std::vector<std::vector<int>> adj(extended_segments.size());
    std::vector<RItem> hits; hits.reserve(64);

    for (int i = 0; i < extended_segments.size(); ++i) {
        // expand query box by r
        BBox q = boxes[i];
        q.min_corner().x(q.min_corner().x() - distance);
        q.min_corner().y(q.min_corner().y() - distance);
        q.max_corner().x(q.max_corner().x() + distance);
        q.max_corner().y(q.max_corner().y() + distance);

        hits.clear();
        rtree.query(bgi::intersects(q), std::back_inserter(hits));

        for (auto const& it : hits) {
            int j = it.second;
            if (j <= i) continue; // add once
            if (CGAL::squared_distance(extended_segments[i].seg, extended_segments[j].seg) > sqd_distance) continue;

            double d = std::fabs(ang[i] - ang[j]);
            if (d > M_PI/2) d = M_PI - d; // fold (treat opposite directions as same)
            if (d > theta_max) continue;

            adj[i].push_back(j);
            adj[j].push_back(i);
        }
    }

    auto groups = group_by_max_degree(adj);
    for (auto &unmodified_id : groups[0]) {
        result[0].push_back(extended_segments[unmodified_id]);
    }

    auto *skip = &groups[0];
    for (auto &id : groups) {
        if (&id==skip)continue;
        std::vector<Segment_w_info> grp;
        for (int i : id) grp.push_back(extended_segments[i]);
        if (!grp.empty()) {
            result.push_back(grp);
        }
    }
    return result;

}
// OUTPUT: groups[0] = all nodes that ended up with degree 0 (never grouped with others)
// groups[1..] = groups formed by repeatedly taking the node of max degree

std::vector<std::vector<int>> group_by_max_degree(const std::vector<std::vector<int>>& adj) {
  const int n = adj.size();
  std::vector<int> deg(n);
  std::vector<char> alive(n, 1);

  // Initialize degrees and find maximum degree
  int Delta = 0;
  for (int u = 0; u < n; ++u) {
    deg[u] = (int)adj[u].size();
    Delta = std::max(Delta, deg[u]);
  }

  // Bucket-queue: bucket[d] holds nodes with degree d (lazy updates allowed)
  std::vector<std::vector<int>> bucket(Delta + 1);
  for (int u = 0; u < n; ++u) bucket[deg[u]].push_back(u);

  // Helper: pop a valid node from the current max-degree bucket
  auto pop_max = [&](int& cur)->int {
    while (cur >= 0) {
      auto &b = bucket[cur];
      while (!b.empty()) {
        int u = b.back(); b.pop_back();
        if (alive[u] && deg[u] == cur) return u; // found a live, up-to-date node
        // otherwise ignore stale entries
      }
      --cur; // move down to next lower degree
    }
    return -1;
  };

  // Result: groups[0] is reserved for all final degree-0 nodes
  std::vector<std::vector<int>> groups;
  groups.emplace_back(); // groups[0] = reserved

  int cur = Delta;
  int alive_cnt = n;

  // Main loop: keep forming groups until the max degree is 0
  while (true) {
    int u = pop_max(cur);
    if (u == -1 || deg[u] == 0) break; // no node with degree > 0 left

    // Build group: u + all currently alive neighbors
    std::vector<int> grp;
    grp.push_back(u);
    for (int v : adj[u]) if (alive[v]) grp.push_back(v);

    // Kill all nodes in this group
    for (int x : grp) {
      if (!alive[x]) continue;
      alive[x] = 0;
      --alive_cnt;
    }

    // Update degrees of remaining alive neighbors
    for (int x : grp) {
      for (int y : adj[x]) if (alive[y]) {
        int d = --deg[y];
        bucket[d].push_back(y); // lazy insert with updated degree
      }
    }

    // Append this group to the solution
    groups.push_back(std::move(grp));
  }

  // Collect all remaining alive nodes (degree 0) into groups[0]
  for (int u = 0; u < n; ++u) if (alive[u]) {
    groups[0].push_back(u);
    alive[u] = 0;
    --alive_cnt;
  }

  return groups;
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





