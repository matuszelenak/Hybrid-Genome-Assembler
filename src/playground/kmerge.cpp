#include <queue>
#include <vector>
#include <iostream>


template<typename T>
void merge_n_vectors(std::vector<std::vector<T> *> &arrays, std::vector<T> &result) {
    std::priority_queue<std::pair<T, int>, std::vector<std::pair<T, int>>, std::greater<std::pair<T, int>>> q;

    int non_exhausted = 0;
    std::vector<int> next_index(arrays.size(), 1);
    for (int i = 0; i < arrays.size(); i++) {
        if (!arrays[i]->empty()) {
            q.push({(*arrays[i])[0], i});
            non_exhausted++;
        }
    }
    if (q.empty()) return;

    auto top_pair = q.top();
    T arr_val = top_pair.first;
    int arr_index = top_pair.second;
    q.pop();

    T previous_value = arr_val;
    result.push_back(previous_value);

    if (next_index[arr_index] == arrays[arr_index]->size()) {
        non_exhausted--;
    } else {
        q.push({(*arrays[arr_index])[next_index[arr_index]], arr_index});
        next_index[arr_index]++;
    }

    while (non_exhausted > 0) {
        top_pair = q.top();
        arr_val = top_pair.first;
        arr_index = top_pair.second;
        q.pop();

        if (arr_val != previous_value) {
            result.push_back(arr_val);
            previous_value = arr_val;
        }

        if (next_index[arr_index] == arrays[arr_index]->size()) {
            non_exhausted--;
        } else {
            q.push({(*arrays[arr_index])[next_index[arr_index]], arr_index});
            next_index[arr_index]++;
        }
    }
}


int main() {
    std::vector<int> a, b, c;
    a = {1, 2, 4, 6, 8, 12, 43};
    b = {2, 3, 5, 5, 7, 12, 53};
    c = {4, 5, 12, 14, 25, 55, 123};

    std::vector<std::vector<int> *> arrays = {&a, &b, &c};
    std::vector<int> merged;
    merge_n_vectors(arrays, merged);

    for (int i = 0; i < merged.size(); i++) {
        std::cout << merged[i] << " ";
    }
}