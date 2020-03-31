#include <string>
#include <fmt/format.h>
#include <iostream>
#include <boost/thread.hpp>

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

void show_progress(uint64_t curr, uint64_t total, const std::string& msg);
int run_command_with_input(const char *command, const std::string& in);


template<typename F, typename T>
auto timeMeasure(const F& func, const T& obj, const std::string& label) {
    return [func, obj, label](auto&&... args) {
        boost::posix_time::ptime ts_start = boost::posix_time::microsec_clock::local_time();
        auto result = (obj->*func)(std::forward<decltype(args)>(args)...);
        std::cout << fmt::format("Execution of {} took {}ms\n", label, (boost::posix_time::microsec_clock::local_time() - ts_start).total_milliseconds());
        return result;
    };
}

template<typename F>
auto timeMeasure(const F& func, const std::string& label) {
    return [func, label](auto&&... args) {
        boost::posix_time::ptime ts_start = boost::posix_time::microsec_clock::local_time();
        auto result = func(std::forward<decltype(args)>(args)...);
        std::cout << fmt::format("Execution of {} took {}ms\n", label, (boost::posix_time::microsec_clock::local_time() - ts_start).total_milliseconds());
        return result;
    };
}

#endif //SRC_UTILS_H
