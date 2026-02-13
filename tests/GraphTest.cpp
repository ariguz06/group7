//
// Created by Patrick Martin on 1/17/26.
//

#include <gtest/gtest.h>
#include "../include/graph/Graph.h"

class Graphs : public testing::Test {
    void SetUp() override {
        adj.reserve(6);

        adj[0].push_back(Edge(1, 1));
        adj[0].push_back(Edge(2, 1));
        adj[0].push_back(Edge(3, 1));
        adj[0].push_back(Edge(4, 1));

        adj[1].push_back(Edge(0, 1));
        adj[1].push_back(Edge(2, 1));

        adj[2].push_back(Edge(0, 1));
        adj[2].push_back(Edge(1, 1));

        adj[3].push_back(Edge(0, 1));

        adj[4].push_back(Edge(0, 1));

        graph = Graph(adj);
    }

protected:
    Graph::AdjMap adj;
    Graph* graph = Graph(adj);
};

TEST_F(Graphs, EdgeExists) {
    graph->populate_buckets();
    EXPECT_TRUE(graph->edge_exists(0, 1));
    EXPECT_TRUE(graph->edge_exists(2, 1));
    EXPECT_TRUE(graph->edge_exists(2, 0));
    EXPECT_TRUE(graph->edge_exists(3, 0));
    EXPECT_TRUE(graph->edge_exists(4, 0));

    EXPECT_FALSE(graph->edge_exists(3, 2));
    EXPECT_FALSE(graph->edge_exists(4, 3));
    EXPECT_FALSE(graph->edge_exists(4, 1));
}

TEST_F(Graphs, PopMinDegree) {
    auto g = graph;

    g->populate_buckets();
    const unsigned long v = g->pop_min_degree_vertex();

    EXPECT_NE(v, 0);
    EXPECT_NE(v, 2);
    EXPECT_NE(v, 1);

    EXPECT_TRUE(v == 4 || v == 3);
}

TEST_F(Graphs, PopVertexNoEdgeFill) {
    const auto g = graph;
    g->populate_buckets();

    unsigned long v = g->pop_min_degree_vertex();
    g->eliminate_vertex(v);

    EXPECT_TRUE(v == 4 || v == 3);

    if (v == 3) {
        EXPECT_FALSE(g->edge_exists(3, 0));
        EXPECT_FALSE(g->edge_exists(0, 3));

        v = g->pop_min_degree_vertex();

        EXPECT_FALSE(v == 3);

    } else {
        EXPECT_FALSE(g->edge_exists(4, 0));
        EXPECT_FALSE(g->edge_exists(0, 4));

        v = g->pop_min_degree_vertex();

        EXPECT_FALSE(v == 4);
    }
}

TEST_F(Graphs, PopVertexEdgeFill) {
    const auto g = graph;
    g->populate_buckets();

    EXPECT_FALSE(g->edge_exists(1, 4));
    EXPECT_FALSE(g->edge_exists(4, 3));
    EXPECT_FALSE(g->edge_exists(3, 2));

    g->eliminate_vertex(0);

    EXPECT_TRUE(g->edge_exists(1, 4));
    EXPECT_TRUE(g->edge_exists(4, 3));
    EXPECT_TRUE(g->edge_exists(3, 2));
}