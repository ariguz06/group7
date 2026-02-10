#pragma once
#include "graph/Graph.h"
#include "util/Timer.h"

#include <queue>
#include <limits>
#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include <sstream>
#include <ctime>

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
    Graph& graph,                // assumes h2h index is already computed 
    unsigned long n,
    unsigned int samples = 10,
    std::ofstream& file,
    Timer& timer
) {

    std::mt19937 rng(42);
    std::uniform_int_distribution<unsigned long> dist(0, n - 1);

    file << "Format <h2h query time>,<dijkstra query time><newline>\n";

    for (unsigned int i = 0; i < samples; i++) {
        unsigned long u = dist(rng);
        unsigned long v = dist(rng);
        
	if (u == v) {
		samples++;
		continue;
	} 


	timer.reset();
	timer.start();
	unsigned long h2h_d = graph.h2h_query(u, v);
	timer.stop();

	file << timer.elapsed() << ",";

	timer.reset();
	timer.start();
        unsigned long oracle_d = dijkstra_oracle(graph, u, v, n);
	timer.stop();

	file << timer.elapsed() << "\n";


        if (h2h_d != oracle_d) {
            // file << " MISMATCH\n";
            // file << "u=" << u << " v=" << v << "\n";
            // file << "H2H=" << h2h_d
            // file << " Oracle=" << oracle_d << "\n";
            return false;
        }
    }

    file << " H2H verified on " << file << " random pairs\n";
    return true;
}

} // namespace GraphUtil


