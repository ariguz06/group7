//
// Created by Patrick Martin on 12/12/25.
//

#ifndef GROUP7_GRAPH_H
#define GROUP7_GRAPH_H
#include <unordered_map>
#include <set>

#include "Edge.h"

class Graph {
    using AdjMap = std::unordered_map<unsigned long, std::vector<Edge>>;
    using TreeDecompAdj = std::unordered_map<unsigned long, std::vector<unsigned long>>; // matrix of edges in TD. key denotes bag root v, edge to bag root u
    using TreeDecompBags = std::unordered_map<unsigned long, AdjMap>; // key: bag root v, value: X(v) represented as adjacency
public:

    explicit Graph(AdjMap adj);

    static Graph from_mtx(const std::string &path, bool weighted=false, bool directed=false);

    static AdjMap get_star(const AdjMap& h, unsigned long vertex);

    static unsigned long treewidth(TreeDecompBags& bags);

    [[nodiscard]] std::vector<unsigned long> bfs_traversal(unsigned long start);

    std::set<unsigned long> get_neighbors(unsigned long vertex);

    void min_degree_elim(AdjMap& g, unsigned long vertex);

    void cache_attributes();

    // returns adj, bags, root of tree decomposition
    std::tuple<TreeDecompAdj, TreeDecompBags, unsigned long> get_td();

private:
    AdjMap adj;

    std::set<std::tuple<unsigned long, unsigned long>> all_edges;

    TreeDecompAdj td_adj;
    TreeDecompBags td_bags;

    unsigned long num_vertices = 0;
};

#endif //GROUP7_GRAPH_H