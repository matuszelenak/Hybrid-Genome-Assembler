#include <cstdio>
#include <cstdlib>
#include <string>
#include "Utils.h"


int run_command_with_input(const char *command, const std::string& in){
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
