#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <optional>
#include <thread>
#include <mutex>

#include "DNAStructures.h"

#ifndef SRC_SEQUENCERECORDITERATOR_H
#define SRC_SEQUENCERECORDITERATOR_H


struct SeqRecordData{
    std::string header;
    std::string sequence;
    std::vector<Quality> qualities;
};


class SequenceRecordIterator {
private:
    std::vector<std::string> paths;
    std::mutex _read_mutex;

    std::ifstream current_file;
    int current_file_index;
    uint64_t current_file_size = 0;
    bool exhausted = false;


    SeqRecordData read_fastq_record(std::string &header);
    SeqRecordData read_fasta_record(std::string &header);
    SeqRecordData (SequenceRecordIterator::*current_record_method)(std::string &header);
public:
    explicit SequenceRecordIterator(std::vector<std::string> &reads_paths);
    ~SequenceRecordIterator();

    std::optional<GenomeReadData> get_next_record();

    std::string get_next_line();

    bool load_file_at_position(int pos);

    std::vector<uint64_t >total_read_bases;
};


#endif //SRC_SEQUENCERECORDITERATOR_H
