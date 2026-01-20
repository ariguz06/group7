//
// Created by Patrick Martin on 12/12/25.
//

#include <utility>
#include <fstream>
#include <sstream>
#include <ranges>
#include <list>
#include <unordered_set>
#include <queue>
#include <unordered_map>
#include <vector>

#include "../../include/graph/Graph.h"

#include <iostream>

Graph::Graph(AdjMap adj) : adj(std::move(adj)) {
    buckets.resize(2000);
    degrees.resize(this->adj.size());
    bucket_position.resize(this->adj.size());
}

Graph Graph::from_mtx(const std::string &path, bool weighted, bool directed) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open file " + path);
    }

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '%') {
            break;
        }
    }

    if (line.empty()) {
        throw std::runtime_error("No matrix header found in file " + path);
    }

    std::istringstream header(line);

    int n, m, nnz;
    header >> n >> m >> nnz;

    AdjMap adj;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '%') continue;
        std::istringstream iss(line);

        unsigned long u, v;
        unsigned long w = 1;
        iss >> u >> v;
        if (weighted)
            iss >> w;

        u--; v--;

        adj[u].push_back({v, w});
        if (!directed && u != v) {
            adj[v].push_back({u, w});
        }
    }

    auto g = Graph(adj);
    g.populate_buckets();

    return g;
}

std::vector<unsigned long> Graph::bfs_traversal(const unsigned long start) {
    std::vector<unsigned long> order;
    std::vector visited(adj.size(), false);

    std::queue<unsigned long> q;
    q.push(start);
    visited[start] = true;

    while (!q.empty()) {
        unsigned long const u = q.front();
        q.pop();
        order.push_back(u);
        for (const auto&[to, w] : adj.at(u)) {
            if (!visited[to]) {
                visited[to] = true;
                q.push(to);
            }
        }

        num_vertices += 1;
    }

    return order;
}

std::vector<unsigned long> Graph::get_neighbors(const unsigned long vertex) {
    std::vector<unsigned long> neighbors;
    for (const auto &[v, w] : adj.at(vertex)) {
        neighbors.push_back(v);
    }

    return neighbors;
}

std::vector<unsigned long> Graph::get_star(const unsigned long vertex) {
    auto star = std::vector<unsigned long>();
    star.push_back(vertex);

    for (const auto &[to, w]: adj.at(vertex)) {
        star.push_back(to);
    }

    return star;
}

void Graph::populate_buckets() {
    for (const auto &[u, edges] : adj) {
        const unsigned long deg = edges.size();

        buckets[deg].push_back(u);
        bucket_position[u] = std::prev(buckets[deg].end());
        degrees[u] = deg;
    }
}

void Graph::eliminate_vertex(const unsigned long v) {
    const auto neighbors = get_neighbors(v);

    // fills in edges w/ updated weights
    for (size_t i = 0; i < neighbors.size(); i++) {
        for (size_t j = i + 1; j < neighbors.size(); j++) {
            unsigned long u = neighbors[i];
            unsigned long w = neighbors[j];

            const unsigned long uvw_weight = get_edge_weight(u, v) + get_edge_weight(v, w);

            if (!edge_exists(u, w)) {
                adj[u].push_back({w, uvw_weight});
                adj[w].push_back({u, uvw_weight});
            } else if (uvw_weight < get_edge_weight(u, w)) {
                adj[u].erase(std::ranges::find_if(adj[u], [&](const Edge& e) {return e.to == w;}));
                adj[w].erase(std::ranges::find_if(adj[w], [&](const Edge& e) {return e.to == u;}));

                adj[u].push_back({w, uvw_weight});
                adj[w].push_back({u, uvw_weight});
            }
        }
    }

    // removes any outward edges from neighbors to vertex
    for (const auto& neighbor : neighbors) {
        adj.at(neighbor).erase(std::ranges::find_if(adj.at(neighbor), [&](const Edge& e) {return e.to == v;}));
    }

    // erases vertex
    adj.erase(v);
    num_vertices -= 1;

    // updates buckets after edge correction
    for (unsigned long neighbor : neighbors) {

        auto& to_erase = adj.at(neighbor);

        const unsigned long d1 = degrees[neighbor];
        const unsigned long d2 = to_erase.size();

        buckets[d1].erase(bucket_position[neighbor]);
        buckets[d2].push_front(neighbor);
        bucket_position[neighbor] = buckets[d2].begin();
        degrees[neighbor] = d2;
    }
}

unsigned long Graph::pop_min_degree_vertex() {
    unsigned long min_bucket = 0;

    while (buckets[min_bucket].empty()) {
        min_bucket++;
    }

    const unsigned long v = buckets[min_bucket].front();
    buckets[min_bucket].pop_front();

    return v;
}

bool Graph::edge_exists(const unsigned long u, const unsigned long v) {
    if (!adj.contains(u)) return false;
    if (!adj.contains(v)) return false;
    return std::ranges::any_of(adj[u].begin(), adj[u].end(), [&](const Edge e) {return e.to == v;});
}

unsigned long Graph::get_edge_weight(const unsigned long u, const unsigned long v) {
    if (u == v) return 0;

    for (const auto& [to, w] : adj[u]) {
        if (to == v) return w;
    }

   return -1;
}

std::tuple<Graph::TreeDecompAdj, Graph::TreeDecompBags, unsigned long> Graph::get_td() {
    auto h = Graph(adj);
    h.num_vertices = 2642;
    h.populate_buckets();

    std::vector<unsigned long> ordering(adj.size());

    td_bags.clear();
    td_adj.clear();

    for (int i = 0; i < adj.size(); i++) {
        unsigned long v = h.pop_min_degree_vertex();

        td_bags[v] = h.get_star(v);
        h.eliminate_vertex(v);

        ordering[v] = i;
    }

    unsigned long root = std::numeric_limits<unsigned long>::max();

    for (unsigned long v = 0; v < ordering.size(); ++v) {
        const auto& bag = td_bags.at(v);

        if (bag.size() == 1) {
            root = v;
            continue;
        }

        unsigned long min_u = v;
        unsigned long best = std::numeric_limits<unsigned long>::max();

        for (unsigned long u : bag) {
            if (u != v && ordering[u] < best) {
                best = ordering[u];
                min_u = u;
            }
        }

        td_adj[v].push_back(min_u);
        td_adj[min_u].push_back(v);
    }

    for (auto &bag: td_bags | std::views::values) {
        std::ranges::sort(bag, [&](const unsigned long a, const unsigned long b) {return ordering[a] > ordering[b];});
    }

    for (auto &[v, neighbors] : td_bags) {
        for (const auto& neighbor : neighbors) {
            td_weights[v].push_back(get_edge_weight(v, neighbor));
        }
    }

    return {td_adj, td_bags, root};
}

unsigned long Graph::treewidth(TreeDecompBags& bags) {
    unsigned long tw = 0;

    for (auto &adj: bags | std::views::values) {
        tw = std::max(tw, adj.size());
    }

    return tw - 1;
}
