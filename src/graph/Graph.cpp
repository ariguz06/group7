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
#include <cassert>
#include <utility>
#include <algorithm>

#include "graph/Graph.h"
#include "util/Oracle.h"

#include <iostream>
#include <random>

Graph::Graph(AdjMap adj, bool populate_buckets) : adj(std::move(adj)) {
    if(populate_buckets) {
        buckets.resize(100000);

        degrees.resize(this->adj.size());
        bucket_position.resize(this->adj.size());
    }
}

// This method generated with Claude Sonnet 4.5
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

    auto g = Graph(adj, false);

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

std::vector<unsigned long> Graph::get_neighbors(const unsigned long vertex) const {
    std::vector<unsigned long> neighbors;
    for (const auto &[v, w] : adj.at(vertex)) {
        neighbors.push_back(v);
    }

    return neighbors;
}

std::vector<unsigned long> Graph::get_star(const unsigned long vertex) const {
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

        for (const auto&[to, w] : edges) {
            add_edge_cache(u, to);
        }
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

            auto& adj_w = adj.at(w);
            auto& adj_u = adj.at(u);

            if (!edge_exists(u, w)) {
                adj_u.push_back({w, uvw_weight});
                adj_w.push_back({u, uvw_weight});

                add_edge_cache(u, w);
                add_edge_cache(w, u);

            } else if (uvw_weight < get_edge_weight(u, w)) {
                // FIXME ?
                adj_u.erase(std::ranges::find_if((adj)[u], [&](const Edge& e) {return e.to == w;}));
                adj_w.erase(std::ranges::find_if(adj_w, [&](const Edge& e) {return e.to == u;}));

                adj_u.push_back({w, uvw_weight});
                adj_w.push_back({u, uvw_weight});
            }
        }
    }

    // add edge weights to td_weights
    for (const auto& neighbor: neighbors) {
        const auto w = get_edge_weight(neighbor, v);
        td_bag_edges[v].push_back({neighbor, w});
    }
    td_bag_edges[v].push_back({v, 0});


    // removes any outward edges from neighbors to vertex
    for (const auto& neighbor : neighbors) {
        remove_edge_cache(neighbor, v);
        remove_edge_cache(v, neighbor);
        adj.at(neighbor).erase(std::ranges::find_if(adj.at(neighbor), [&](const Edge& e) {return e.to == v;}));
    }

    // erases vertex
    adj.erase(v);
    num_vertices -= 1;

    // updates buckets after edge fill-in
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

// obtains vertex with minimum degree and pops it from its corresponding bucket
unsigned long Graph::pop_min_degree_vertex() {
    unsigned long min_bucket = 0;

    while (buckets[min_bucket].empty()) {
        min_bucket++;
    }

    const unsigned long v = buckets[min_bucket].front();
    buckets[min_bucket].pop_front();

    return v;
}

bool Graph::edge_exists(const unsigned long u, const unsigned long v) const {
    const uint64_t edge = (static_cast<uint64_t>(u) << 32) | static_cast<uint64_t>(v);
    return edge_set.contains(edge);
}

unsigned long Graph::get_edge_weight(const unsigned long u, const unsigned long v) const {
    if (u == v) return 0;

    for (const Edge& e : adj.at(u)) {
        if (e.to == v) return e.w;
    }

    return ULONG_MAX;
}

void Graph::add_edge_cache(const unsigned long u, const unsigned long v) {
    // first 32 bits encode u (first shifted left 32 bits), last 32 bits encode v
    const uint64_t edge_cache = (static_cast<uint64_t>(u) << 32) | static_cast<uint64_t>(v);

    edge_set.insert(edge_cache);
}

void Graph::remove_edge_cache(const unsigned long u, const unsigned long v) {
    const uint64_t edge_cache = (static_cast<uint64_t>(u) << 32) | static_cast<uint64_t>(v);

    edge_set.erase(edge_cache);
}

std::vector<unsigned long> Graph::get_random_ordering() const {
    std::vector<unsigned long> ordering(adj.size());
    for (unsigned long u = 0; u < adj.size(); u++) {
        ordering[u] = u;
    }

    std::mt19937_64 rng(std::random_device{}());
    std::ranges::shuffle(ordering, rng);

    return ordering;
}

std::tuple<Graph::TreeDecompAdj, Graph::TreeDecompBags, unsigned long> Graph::get_td() {
    Graph h = Graph(adj, true);
    h.num_vertices = adj.size();
    h.populate_buckets();

    std::vector<unsigned long> ordering(adj.size());

    td_bags.clear();
    td_adj.clear();

    // const auto rand_ordering = get_random_ordering();

    for (int i = 0; i < adj.size(); i++) {
        unsigned long v = h.pop_min_degree_vertex(); // can be substituted with other heuristic
        // unsigned long v = rand_ordering[i];

        td_bags[v] = h.get_star(v);
        h.eliminate_vertex(v);

        ordering[v] = i;

        if (i % static_cast<int>(num_vertices / 10) == 0) {
            std::cout << "Eliminated vertex " << i << " " << 10 * i / static_cast<int>(num_vertices / 10) << "%" << std::endl;
        }
    }

    for (unsigned long v : adj | std::views::keys) {
        const auto& bag = td_bags.at(v);

        unsigned long min_u = v;
        unsigned long best = ULONG_MAX;

        for (const unsigned long u : bag) {
            if (u != v && ordering[u] < best) {
                best = ordering[u];
                min_u = u;
            }
        }

        if (best == ULONG_MAX) {
            td_root = v;
            parent_map[v] = v;
            continue;
        }

        td_adj[v].push_back(min_u);
        td_adj[min_u].push_back(v);

        parent_map[v] = min_u;
    }

    for (auto& [v, bag] : td_bags) {

        std::ranges::sort(bag, [&](const unsigned long a, const unsigned long b) {return ordering[a] > ordering[b];});

        if (v == td_root) {
            td_weights[v].push_back(0);
            continue;
        }

        const auto& edges = h.td_bag_edges.at(v);

        for (const unsigned long u : bag) {
            
            for (const auto&[to, w] : edges) {
                if (to == u && u != v) {
                    td_weights[v].push_back(w);
                }
            }
        }

        td_weights[v].push_back(0);
    }

    return {td_adj, td_bags, td_root};
}

std::vector<unsigned long> Graph::get_top_down_ordering() const {
    std::vector<unsigned long> ordering;

    std::queue<unsigned long> q;
    std::unordered_set<unsigned long> visited;

    q.push(td_root);
    visited.insert(td_root);

    while (!q.empty()) {
        unsigned long u = q.front();
        q.pop();
        ordering.push_back(u);

        for (const unsigned long neighbor : td_adj.at(u)) {
            if (!visited.contains(neighbor)) {
                q.push(neighbor);
                visited.insert(neighbor);
            }
        }
    }

    return ordering;
}

std::vector<unsigned long> Graph::get_bag_path(const unsigned long v) const {
    std::vector<unsigned long> path;

    if (v == td_root) {
        path.push_back(v);
        return path;
    }

    unsigned long current = v;
    while (current != td_root) {
        path.push_back(current);
        
        unsigned long parent = parent_map.at(current);
        current = parent;
    }
    
    path.push_back(td_root);
    std::ranges::reverse(path);
    
    return path;
}

unsigned long Graph::lca(const unsigned long u, const unsigned long v) {
    const auto& u_anc = anc_map[u];
    const auto& v_anc = anc_map[v];

    const auto min_len = std::min(u_anc.size(), v_anc.size());
    unsigned long lca = ULONG_MAX;

    for (int i = 0; i < min_len; i++) {
        if (u_anc[i] == v_anc[i]) {
            lca = u_anc[i];
        }
        else {
            break;
        }
    }

    if (lca == ULONG_MAX) {
        throw std::invalid_argument("Graph::lca() failed");
    }

    return lca;
}

unsigned long Graph::h2h_query(const unsigned long u, const unsigned long v) {

    const auto x = lca(u, v);

    unsigned long d = 1e9;
    const auto& dis_map = std::get<1>(*h2h);
    const auto& pos_map = std::get<0>(*h2h);

    for (const unsigned long i : pos_map.at(x)) {
        d = std::min(d, dis_map.at(u)[i] + dis_map.at(v)[i]);
    }

    return d;
}

unsigned long Graph::index_of(const std::vector<unsigned long>& b, const unsigned long v) {
    for (unsigned long i = 0; i < b.size(); ++i) {
        if (b[i] == v) {return i;}
    }

    throw std::out_of_range("Vertex not found");
}

std::tuple<Graph::Pos, Graph::Dis> Graph::get_h2h() {
    Pos pos;
    Dis dis;

    const std::vector<unsigned long> ordering = get_top_down_ordering();

    for (const unsigned long bag : ordering) {
        anc_map[bag] = get_bag_path(bag);
    }

    for (const unsigned long v_bag : ordering) {

        auto& anc = anc_map.at(v_bag);

        for (const unsigned long bag_vertex : td_bags.at(v_bag)) {
            auto bag_pos_i = index_of(anc, bag_vertex);
            pos[v_bag].push_back(bag_pos_i);
        }
    }

    for (std::vector<unsigned long>& pos_arr : pos | std::views::values) {
        std::ranges::sort(pos_arr);
    }

    for (const unsigned long v_bag : ordering) {
        auto& anc = anc_map.at(v_bag);
        auto& bag = td_bags.at(v_bag);    

        for (unsigned long i = 0; i < anc.size()-1; i++) {

            dis[v_bag].push_back(1e9);

            for (unsigned long j = 0; j < bag.size()-1; j++) {
                unsigned long d;
                unsigned long xj_vertex = bag.at(j);

                if (pos[v_bag][j] > i) {
                    d = dis[xj_vertex][i];
                } else {
                    unsigned long v_bag_anc_i = anc_map.at(v_bag)[i];
                    const unsigned long dis_index = pos[v_bag][j];

                    d = dis[v_bag_anc_i][dis_index];
                }

                dis[v_bag][i] = std::min(dis[v_bag][i], td_weights.at(v_bag)[j] + d);
            }
        }

        dis[v_bag].push_back(0);
    }

    h2h = std::make_unique<std::tuple<Pos, Dis>>(std::make_tuple(pos, dis));
    return *h2h;
}


unsigned long Graph::treewidth(TreeDecompBags& bags) {
    unsigned long tw = 0;

    for (auto &bag: bags | std::views::values) {
        if (bag.size() > tw) {
            tw = bag.size();
        }
    }

    return tw - 1;
}
