#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <fmt/format.h>

#include "Utils.h"


void progress_bar(uint64_t curr, uint64_t total, const std::string& msg) {
    double progress_percent = (double) curr / (double) total;

    std::cout << std::fixed << std::setprecision(2)
              << fmt::format("\r{}: [{}{}] {}%", msg, std::string(progress_percent * 100, '#'), std::string(100 - progress_percent * 100, ' '), (int) (100 * progress_percent));

    if (progress_percent == 1)
        std::cout << std::endl;
    else
        std::cout.flush();
}


int run_command_with_input(const char *command, const std::string &in) {
    FILE *stream;

    stream = popen(command, "w");
    if (!stream) {
        fprintf(stderr,
                "incorrect parameters or too many files.\n");
        return EXIT_FAILURE;
    }

    fprintf(stream, "%s", in.c_str());
    if (ferror(stream)) {
        fprintf(stderr, "Output to stream failed.\n");
        exit(EXIT_FAILURE);
    }

    if (pclose(stream) != 0) {
        fprintf(stderr,
                "Could not run more or other error.\n");
    }
    return EXIT_SUCCESS;
}

std::string capture_output_of_command(const char *command){
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(command, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

std::vector<std::string> split_string(std::string &s, std::string delim){
    std::vector<std::string> result;
    auto start = 0U;
    auto end = s.find(delim);
    while (end != std::string::npos)
    {
        result.push_back(s.substr(start, end - start));
        start = end + delim.length();
        end = s.find(delim, start);
    }
    return result;
}