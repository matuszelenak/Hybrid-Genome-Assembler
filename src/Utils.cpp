#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <fmt/format.h>
#include "Utils.h"


void show_progress(int curr, int total, std::string msg) {
    float progress_percent = (float) curr / (float) total;

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
