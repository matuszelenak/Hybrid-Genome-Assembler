#include <string>
#include <vector>
#include <fmt/format.h>
#include <iostream>
#include <boost/thread.hpp>
#include <thread>
#include <queue>

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

class ThreadRunner {
public:

    unsigned int num_of_threads;
    std::vector<std::thread *> threads;

    ThreadRunner() {
        num_of_threads = std::thread::hardware_concurrency();
    };

    ~ThreadRunner() {
        for (auto thread: threads) free(thread);
    }

    explicit ThreadRunner(auto &&... args) : ThreadRunner() {
        for (int i = 0; i < num_of_threads; i++) {
            add_thread(std::forward<decltype(args)>(args)...);
        }
        run();
    }

    bool add_thread(auto &&... args) {
        if (threads.size() == num_of_threads) return false;
        threads.push_back(new std::thread(std::forward<decltype(args)>(args)...));
        return true;
    };

    void run() { for (auto thread : threads) thread->join(); };
};

template<typename T>
void merge_n_vectors(std::vector<std::vector<T> *> &arrays, std::vector<T> &result, bool unique) {
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
}

#endif //SRC_UTILS_H
