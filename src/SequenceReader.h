//
// Created by whiskas on 04/02/19.
//

#include <string>
#include <fstream>
#include <stdexcept>
#include <optional>

#ifndef SRC_SEQUENCEREADER_H
#define SRC_SEQUENCEREADER_H


class SequenceReader {
private:
    std::ifstream input_file;

    std::optional<std::string> read_fasta_sequence();
    std::optional<std::string> read_fastq_sequence();
public:
    explicit SequenceReader(const std::string &path);
    std::optional<std::string> (SequenceReader::*read_sequence_line)();
    std::optional<std::string> get_next_record();
};


#endif //SRC_SEQUENCEREADER_H
