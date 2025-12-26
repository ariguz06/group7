#include <iostream>

#include "include/graph/Graph.h"
#include "include/util/Timer.h"

// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
int main() {

    Graph graph = Graph::from_mtx("mtx/road-usroads-48.mtx", false, false);
    auto t = Timer();

    t.start();
    const std::vector<unsigned long> order = graph.bfs_traversal(7);
    t.stop();

    std::cout << "BFS time elapsed: " << t.elapsed() << " with size " << order.size() << std::endl;
    t.reset();

    t.start();
    auto td_metrics = graph.get_td();
    t.stop();

    std::cout << "Tree decomposition time elapsed: " << t.elapsed() << std::endl;
    std::cout << "Treewidth: " << Graph::treewidth(std::get<1>(td_metrics));

    return 0;
}