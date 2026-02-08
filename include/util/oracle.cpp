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

        if (u == target) return d;
        if (d > dist[u]) continue;

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
    Graph& graph,                // ❗ NOT const
    unsigned long n,
    unsigned int samples = 100
) {
    // Build H2H once
    graph.get_h2h();

    std::mt19937 rng(42);
    std::uniform_int_distribution<unsigned long> dist(0, n - 1);

    for (unsigned int i = 0; i < samples; i++) {
        unsigned long u = dist(rng);
        unsigned long v = dist(rng);
        if (u == v) continue;

        unsigned long h2h_d = graph.h2h_query(u, v);
        unsigned long oracle_d = dijkstra_oracle(graph, u, v, n);

        if (h2h_d != oracle_d) {
            std::cout << " MISMATCH\n";
            std::cout << "u=" << u << " v=" << v << "\n";
            std::cout << "H2H=" << h2h_d
                      << " Oracle=" << oracle_d << "\n";
            return false;
        }
    }

    std::cout << " H2H verified on " << samples << " random pairs\n";
    return true;
}

} // namespace GraphUtil

int main() {
    Graph graph = Graph::from_mtx(
        "road-minnesota.mtx",
        true,   // weighted
        false   // undirected
    );

    auto [td_adj, td_bags, root] = graph.get_td();

    std::cout << "Treewidth: "
              << Graph::treewidth(td_bags)
              << std::endl;

    constexpr unsigned long N = 2642; // road-minnesota size

    GraphUtil::verify_h2h(graph, N, 50);

    return 0;
}
