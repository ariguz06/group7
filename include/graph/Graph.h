//
// Created by Patrick Martin on 12/12/25.
//

#ifndef GROUP7_GRAPH_H
#define GROUP7_GRAPH_H
#include <list>
#include <unordered_map>
#include <unordered_set>

#include "Edge.h"

class Graph {

public:

    using AdjMap = std::unordered_map<unsigned long, std::vector<Edge>>;

    using TreeDecompAdj = std::unordered_map<unsigned long, std::vector<unsigned long>>; // matrix of edges in TD. key denotes bag root v, edge to bag root u
    using TreeDecompBags = TreeDecompAdj; // key: bag root v, value: all vertices in X(v) but without weights

    using Pos = std::unordered_map<unsigned long, std::vector<unsigned long>>;
    using Dis = std::unordered_map<unsigned long, std::vector<unsigned long>>;

    // key v denotes X(v), value psi_i, where 0 <= i < |X(v)|  denotes weight of each vertex in X(v) to v
    // populated after sorting bags by ordering
    using TreeDecompWeights = std::unordered_map<unsigned long, std::vector<unsigned long>>;

    Graph() = default;

    explicit Graph(AdjMap adj);

    static Graph from_mtx(const std::string &path, bool weighted=false, bool directed=false);

    [[nodiscard]] std::vector<unsigned long> bfs_traversal(unsigned long start);

    // returns adj, bags, root of tree decomposition
    std::tuple<TreeDecompAdj, TreeDecompBags, unsigned long> get_td();

    [[nodiscard]] std::tuple<Pos, Dis> get_h2h();
    [[nodiscard]] std::vector<unsigned long> get_top_down_ordering() const;
    [[nodiscard]] std::vector<unsigned long> get_bag_path(unsigned long v) const;
    [[nodiscard]] unsigned long h2h_query(unsigned long u, unsigned long v);

    void populate_buckets();

    [[nodiscard]] bool edge_exists(unsigned long u, unsigned long v) const;
    [[nodiscard]] unsigned long pop_min_degree_vertex();
    [[nodiscard]] std::vector<unsigned long> get_star(unsigned long vertex) const;
    void eliminate_vertex(unsigned long v);
    static unsigned long treewidth(TreeDecompBags& bags);

    [[nodiscard]] std::vector<unsigned long> get_neighbors(unsigned long vertex) const;
    [[nodiscard]] unsigned long get_edge_weight(unsigned long u, unsigned long v) const;

    void add_edge_cache(unsigned long u, unsigned long v);
    void remove_edge_cache(unsigned long u, unsigned long v);

    [[nodiscard]] std::vector<unsigned long> get_random_ordering() const;

private:
    AdjMap* adj{};
    unsigned long num_vertices = 0;

    TreeDecompAdj td_adj;
    TreeDecompBags td_bags;
    TreeDecompWeights td_weights;

    std::tuple<Pos, Dis> h2h;

    std::unordered_map<unsigned long, std::vector<unsigned long>> anc_map;

    unsigned long td_root = std::numeric_limits<unsigned long>::max();

    std::unordered_set<uint64_t> edge_set;

    std::vector<std::list<unsigned long>> buckets;
    std::vector<unsigned long> degrees;
    std::vector<std::list<unsigned long>::iterator> bucket_position;

    static unsigned long index_of(const std::vector<unsigned long>&, unsigned long v);

    unsigned long lca(unsigned long u, unsigned long v);
};

#endif //GROUP7_GRAPH_H