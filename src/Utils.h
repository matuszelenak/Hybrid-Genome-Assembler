#include <string>

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

void show_progress(uint64_t curr, uint64_t total, const std::string& msg);
int run_command_with_input(const char *command, const std::string& in);

#endif //SRC_UTILS_H
