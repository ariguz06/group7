#include <filesystem>
#include <iostream>
#include <chrono>
#include <sstream>
#include <fstream>
#include <ctime>

#include "graph/Graph.h"
#include "util/Timer.h"
// #include "util/Oracle.h"

#define OUTPUT_DIR "/out"

int main() {

    auto t = Timer();

    std::vector<std::string> file_paths;

    for (const auto& filepath : std::filesystem::directory_iterator("./mtx")) {
        file_paths.push_back(filepath.path());
    }

    auto now = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&time);

    std::ostringstream oss;
    oss << std::setfill('0') 
	    << std::setw(2) << tm.tm_mon + 1
	    << std::setw(2) << tm.tm_mday
	    << std::setw(2) << tm.tm_hour
	    << std::setw(2) << tm.tm_min;

    for (const auto& filepath : file_paths) {
		std::string path = "./out/" + oss.str() + ".txt";

		std::cout << "Creating file at location " << path << std::endl;

		std::ofstream file(path);

		if (!file.is_open()) {
			std::cerr << "Failed to open file" << std::endl;
			return 1;
		}

		t.start();
		Graph graph = Graph::from_mtx(filepath, false, false);
		t.stop();

		file << "Graph construction time elapsed: " << std::to_string(t.elapsed()) << "\n\n";
		t.reset();

		t.start();
		auto td_metrics = graph.get_td();
		t.stop();

		auto& td_bags = std::get<1>(td_metrics);

		file << "TD time elapsed: " << std::to_string(t.elapsed()) << "\n";
		file << "Treewidth: " << std::to_string(t.elapsed()) << "\n\n";
		t.reset();

		t.start();
		graph.get_h2h();
		t.stop();

		file << "H2H construction time elapsed : " << std::to_string(t.elapsed()) << "\n";
		file << "H2H size: " << std::to_string(graph.get_h2h_size()) << "\n";

		Oracle::verify_h2h(graph, graph.num_vertices, 10, file, t);

		file << "End" << "\n";
		file.close();
    }

/* Benchmarking plan to be run on Unity Cluster w/ job script

TD benchmark - compute TDs + width of every graph
H2H construction time + shortest distance query times (random sample)
Dijkstra's shortest distance query times (random sample)

*/

    return 0;
}
