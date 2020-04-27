#include <string>
#include <vector>
#include <fmt/format.h>
#include <iostream>
#include <boost/thread.hpp>
#include <thread>
#include <queue>
#include <mutex>

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

void progress_bar(uint64_t curr, uint64_t total, const std::string &msg);

int run_command_with_input(const char *command, const std::string &in);


template<typename F, typename T>
auto timeMeasureMemberFunc(const F &func, const T &obj, const std::string &label) {
    return [func, obj, label](auto &&... args) {
        boost::posix_time::ptime ts_start = boost::posix_time::microsec_clock::local_time();
        auto result = (obj->*func)(std::forward<decltype(args)>(args)...);
        std::cout << fmt::format("Execution of {} took {}ms\n", label, (boost::posix_time::microsec_clock::local_time() - ts_start).total_milliseconds());
        return result;
    };
}

template<typename F>
auto timeMeasure(const F &func, const std::string &label) {
    return [func, label](auto &&... args) {
        boost::posix_time::ptime ts_start = boost::posix_time::microsec_clock::local_time();
        auto result = func(std::forward<decltype(args)>(args)...);
        std::cout << fmt::format("Execution of {} took {}ms\n", label, (boost::posix_time::microsec_clock::local_time() - ts_start).total_milliseconds());
        return result;
    };
}

void run_in_threads(auto &&... args){
    unsigned int num_of_threads = std::thread::hardware_concurrency();
    std::thread threads[num_of_threads];
    for (unsigned int i = 0; i < num_of_threads; i++){
        threads[i] = std::thread(std::forward<decltype(args)>(args)...);
    }
    for (unsigned int i = 0; i < num_of_threads; i++){
        threads[i].join();
    }
}

template<typename T>
class ConcurrentQueue {
    std::mutex mut;
public:
    std::queue<T> queue;
    ConcurrentQueue() = default;;

    template<typename It>
    explicit ConcurrentQueue(It begin, It end){
        for (It it = begin; it != end; ++it){
            queue.push(*it);
        }
    }

    std::optional<T> pop(){
        std::lock_guard<std::mutex> lock(mut);
        if (!queue.empty()){
            auto result = std::optional<T>{queue.front()};
            queue.pop();
            return result;
        }
        return std::nullopt;
    }
};

template<typename T>
std::vector<T> merge_n_vectors(std::vector<std::vector<T> *> &arrays, bool unique) {
    std::priority_queue<std::pair<T, int>, std::vector<std::pair<T, int>>, std::greater<std::pair<T, int>>> q;

    int non_exhausted = 0;
    std::vector<int> next_index(arrays.size(), 1);
    for (int i = 0; i < arrays.size(); i++) {
        if (!arrays[i]->empty()) {
            q.push({(*arrays[i])[0], i});
            non_exhausted++;
        }
    }
    if (q.empty()) return {};

    auto top_pair = q.top();
    T arr_val = top_pair.first;
    int arr_index = top_pair.second;
    q.pop();

    std::vector<T> result;
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

        if (!unique || arr_val != previous_value) {
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
    return result;
}

template<typename T>
std::vector<T> get_vectors_intersection(std::vector<T> &x, std::vector<T> &y){
    std::vector<T> result;
    int i = 0, j = 0;
    while (i < x.size() && j < y.size()) {
        if (x[i] < y[j])
            i++;
        else if (y[j] < x[i])
            j++;
        else
        {
            result.push_back(x[i]);
            i++;
            j++;
        }
    }
    return result;
}

#endif //SRC_UTILS_H
