//
// Created by whiskas on 04/02/19.
//

#include <string>
#include <fstream>
#include <stdexcept>
#include <optional>
#include <thread>
#include <mutex>

#include "structs.h"

#ifndef SRC_SEQUENCEREADER_H
#define SRC_SEQUENCEREADER_H


class SequenceReader {
private:
    std::ifstream input_file;

    std::optional<GenomeRead> read_fasta_sequence();
    std::optional<GenomeRead> read_fastq_sequence();

    std::mutex _read_mutex;
public:
    explicit SequenceReader(const std::string &path);
    std::optional<GenomeRead> (SequenceReader::*read_sequence_line)();
    std::optional<GenomeRead> get_next_record();
};


#endif //SRC_SEQUENCEREADER_H
