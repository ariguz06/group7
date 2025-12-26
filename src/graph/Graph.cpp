//
// Created by Patrick Martin on 12/12/25.
//

#include <utility>
#include <fstream>
#include <sstream>

#include "../../include/graph/Graph.h"
#include "../../include/util/MinHeap.h"

#include <iostream>
#include <ranges>
#include <unordered_set>

Graph::Graph(AdjMap adj) : adj(std::move(adj)) {}

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

    auto to_return = Graph(adj);
    to_return.cache_attributes();

    return to_return;
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

void Graph::cache_attributes() {
    for (const auto &[u, edges] : adj) {
        for (const auto &[v, w] : edges) {
            auto edge = std::make_tuple(u, v);
            all_edges.insert(edge);
        }
    }
}

std::set<unsigned long> Graph::get_neighbors(const unsigned long vertex) {
    std::set<unsigned long> neighbors;
    for (const auto &[v, w] : adj.at(vertex)) {
        if (!neighbors.contains(v)) {
            neighbors.insert(v);
        }
    }

    return neighbors;
}

Graph::AdjMap Graph::get_star(const AdjMap& h, const unsigned long vertex) {
    if (!h.contains(vertex)) {
        throw std::runtime_error("Vertex " + std::to_string(vertex) + " does not exist in graph");
    }

    auto map = std::unordered_map<unsigned long, std::vector<Edge>>();
    map[vertex] = h.at(vertex);
    map[vertex].push_back({vertex, 0});

    return map;
}

void Graph::min_degree_elim(AdjMap &g, const unsigned long vertex) {
    const std::set<unsigned long> neighbors = get_neighbors(vertex);

    for (const auto u : neighbors) {
        for (const auto w : neighbors) {

            auto weight_u_v_it = std::ranges::find_if(adj.at(u),
                                                      [vertex](const Edge &e) { return e.to == vertex; });
            if (weight_u_v_it == adj.at(u).end()) {
                throw std::runtime_error("Error finding weight from u to vertex");
            }
            const unsigned long weight_u_v = weight_u_v_it->w;


            auto weight_v_w_it = std::ranges::find_if(adj.at(vertex),
                                                      [w](const Edge &e) { return e.to == w; });
            if (weight_v_w_it == adj.at(vertex).end()) {
                throw std::runtime_error("Error finding weight from vertex to w");
            }
            const unsigned long weight_v_w = weight_v_w_it->w;

            auto existing_edge_it = std::ranges::find_if(g[u],
                                                         [w](const Edge &e) { return e.to == w; });

            if (existing_edge_it == g[u].end()) {
                g[u].push_back({w, weight_u_v + weight_v_w});
                g[w].push_back({u, weight_u_v + weight_v_w});
                all_edges.insert(std::make_tuple(u, w));
                all_edges.insert(std::make_tuple(w, u));
            } else {
                if (weight_u_v + weight_v_w < existing_edge_it->w) {
                    existing_edge_it->w = weight_u_v + weight_v_w;

                    auto reverse_edge_it = std::ranges::find_if(g[w],
                                                                [u](const Edge &e) { return e.to == u; });
                    if (reverse_edge_it != g[w].end()) {
                        reverse_edge_it->w = weight_u_v + weight_v_w;
                    }
                }
            }
        }
    }

    g.erase(vertex);
    for (auto neighbor : neighbors) {
        auto &outward_edges = g[neighbor];
        std::erase_if(outward_edges, [vertex](const Edge &e){ return e.to == vertex; });
    }
}

std::tuple<Graph::TreeDecompAdj, Graph::TreeDecompBags, unsigned long> Graph::get_td() {
    auto h = adj;

    std::vector<unsigned long> ordering(adj.size());

    unsigned long root = 0;
    unsigned long step = 0;
    td_bags.clear();
    td_adj.clear();

    MinHeap minHeap;

    for (const auto &[v, edges]: adj) {
        minHeap.insert(v, edges.size());
    }

    for (int i = 0; i < adj.size(); i++) {
        unsigned long v_min_degree = minHeap.extract_min();

        td_bags[v_min_degree] = get_star(h, v_min_degree);

        ordering[v_min_degree] = step++;
        min_degree_elim(h, v_min_degree);

        if (i % 1000 == 0) {
            std::cout << "Finished bag " << i << std::endl;
        }
    }

    std::cout << "Finished ordering" << std::endl;

    for (const auto& [v, bag] : td_bags) {
        if (bag.empty()) continue;

        unsigned long min_ordering_val = std::numeric_limits<unsigned long>::max();
        unsigned long min_ordering_vertex = v;

        for (const auto &vertex_id: bag | std::views::keys) {
            if (vertex_id == v) continue;
            if (ordering[vertex_id] < min_ordering_val) {
                min_ordering_val = ordering[vertex_id];
                min_ordering_vertex = vertex_id;
            }
        }

        if (min_ordering_vertex != v) {
            td_adj[min_ordering_vertex].push_back(v);
            td_adj[v].push_back(min_ordering_vertex);
            root = min_ordering_vertex;
        }
    }

    return {td_adj, td_bags, root};
}

unsigned long Graph::treewidth(TreeDecompBags& bags) {
    unsigned long tw = 0;

    for (auto &adj: bags | std::views::values) {
        for (auto &vertices : adj | std::views::values) {
            tw = std::max(tw, vertices.size());
        }
    }

    return tw - 1;
}


