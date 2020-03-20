#include <string>
#include <fstream>
#include <stdexcept>
#include <optional>
#include <thread>
#include <mutex>

#include "DNAStructures.h"
#include "Utils.h"

#ifndef SRC_READDATALOADER_H
#define SRC_READDATALOADER_H


class ReadDataLoader {
private:
    std::ifstream input_file;
    TwoCategoriesIdentifier categories = TwoCategoriesIdentifier();

    std::optional<GenomeReadData> read_fasta_record();
    std::optional<GenomeReadData> read_fastq_record();

    std::mutex _read_mutex;

    std::string get_header();
public:
    explicit ReadDataLoader(const std::string &path);
    ~ReadDataLoader();

    std::optional<GenomeReadData> (ReadDataLoader::*read_record)();

    std::optional<GenomeReadData> get_next_record();
};


#endif //SRC_READDATALOADER_H
