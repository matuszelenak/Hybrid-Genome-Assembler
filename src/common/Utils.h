#include <string>
#include <fmt/format.h>
#include <iostream>
#include <boost/thread.hpp>
#include <thread>

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

void show_progress(uint64_t curr, uint64_t total, const std::string &msg);

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

    ThreadRunner(auto &&... args) : ThreadRunner() {
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

#endif //SRC_UTILS_H
