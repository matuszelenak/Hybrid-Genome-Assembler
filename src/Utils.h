#include <string>
#include <map>

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

class TwoCategoriesIdentifier {
    std::map<std::string, bool> _category_to_id;
public:
    void add_category(const std::string& category);
    bool id_for_category(const std::string& category);
};

int run_command_with_input(const char *command, const std::string& in);

#endif //SRC_UTILS_H
