#pragma once
#include "../graph/Graph.h"
#include <queue>
#include <limits>
#include <vector>
#include <random>
#include <iostream>

namespace GraphUtil {

inline unsigned long dijkstra_oracle(
    const Graph& graph,
    unsigned long source,
    unsigned long target,
    unsigned long n
) {
    const unsigned long INF = std::numeric_limits<unsigned long>::max();

    std::vector<unsigned long> dist(n, INF);
    dist[source] = 0;

    using P = std::pair<unsigned long, unsigned long>; // (dist, node)
    std::priority_queue<P, std::vector<P>, std::greater<>> pq;
    pq.push({0, source});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (d > dist[u]) continue;
        if (u == target) return dist[target];

        for (unsigned long v : graph.get_neighbors(u)) {
            unsigned long w = graph.get_edge_weight(u, v);
            if (dist[v] > d + w) {
                dist[v] = d + w;
                pq.push({dist[v], v});
            }
        }
    }
    return dist[target];
}

inline bool verify_h2h(
    Graph& graph,                // 
    unsigned long n,
    unsigned int samples,
    std::ofstream& file
) {
    // Build H2H once

    std::mt19937 rng(42);
    std::uniform_int_distribution<unsigned long> dist(0, n - 1);

    bool valid = true;

    for (unsigned int i = 0; i < samples; i++) {
        unsigned long u = dist(rng);
        unsigned long v = dist(rng);
        if (u == v) { samples++; continue; };

        unsigned long h2h_d = graph.h2h_query(u, v);
        unsigned long oracle_d = dijkstra_oracle(graph, u, v, n);

        if (h2h_d != oracle_d) {
            std::cout << " MISMATCH\n";
            std::cout << "u=" << u << " v=" << v << "\n";
            std::cout << "H2H=" << h2h_d
                      << " Oracle=" << oracle_d << "\n";

            valid = false;
        }
    }

    return valid;
}

} // namespace GraphUtil
