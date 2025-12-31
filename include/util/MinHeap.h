//
// Created by Patrick Martin on 12/16/25.
//


// MinHeap.h
#ifndef MINHEAP_H
#define MINHEAP_H

#include <vector>
#include <unordered_map>

class MinHeap {

    std::vector<std::pair<unsigned long, unsigned long>> heap;

    std::unordered_map<unsigned long, size_t> id_to_index;

    void heapify_up(size_t index) {
        while (index > 0 && heap[index].second < heap[(index - 1) / 2].second) {
            swap_elements(index, (index - 1) / 2);
            index = (index - 1) / 2;
        }
    }

    void heapify_down(const size_t index) {
        size_t smallest = index;
        const size_t left = 2 * index + 1;
        const size_t right = 2 * index + 2;

        if (left < heap.size() && heap[left].second < heap[smallest].second) {
            smallest = left;
        }
        if (right < heap.size() && heap[right].second < heap[smallest].second) {
            smallest = right;
        }

        if (smallest != index) {
            swap_elements(index, smallest);
            heapify_down(smallest);
        }
    }

    void swap_elements(size_t i, size_t j) {
        std::swap(heap[i], heap[j]);
        id_to_index[heap[i].first] = i;
        id_to_index[heap[j].first] = j;
    }

public:

    void reserve(const size_t size) {
        heap.reserve(size);
    }

    void insert(unsigned long vertex, unsigned long degree) {
        heap.emplace_back(vertex, degree);
        id_to_index[vertex] = heap.size() - 1;
        heapify_up(heap.size() - 1);
    }

    // Update the degree of a vertex in the heap
    void update(const unsigned long vertex, unsigned long new_degree) {
        if (!id_to_index.contains(vertex)) {
            return; // Vertex not found
        }
        const size_t index = id_to_index[vertex];
        heap[index].second = new_degree;
        heapify_up(index);
        heapify_down(index);
    }

    unsigned long extract_min() {
        if (heap.empty()) {
            throw std::runtime_error("Heap is empty");
        }

        const unsigned long min_vertex = heap[0].first;
        swap_elements(0, heap.size() - 1);
        heap.pop_back();
        id_to_index.erase(min_vertex);

        if (!heap.empty()) {
            heapify_down(0);
        }
        return min_vertex;
    }

    // Check if the heap is empty
    [[nodiscard]] bool empty() const {
        return heap.empty();
    }
};

#endif // MINHEAP_H