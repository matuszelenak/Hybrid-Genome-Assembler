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


struct SeqRecordData {
    std::string header;
    std::string sequence;
    std::vector<Quality> qualities;
};


struct ReadFileMetaData {
    std::string filename;
    uint64_t records;
    uint64_t min_read_length;
    uint64_t max_read_length;
    uint64_t avg_read_length;
    uint64_t total_bases;
};


class SequenceRecordIterator {
private:
    std::vector<std::string> paths;
    std::mutex _read_mutex;

    std::ifstream current_file;
    int current_file_index = 0;
    uint64_t current_file_size = 0;
    bool exhausted = false;


    SeqRecordData read_fastq_record(std::string &header);

    SeqRecordData read_fasta_record(std::string &header);

    SeqRecordData (SequenceRecordIterator::*current_record_method)(std::string &header) = &SequenceRecordIterator::read_fastq_record;

    std::string get_next_line();

    bool load_file_at_position(int pos);

    void get_meta_data();
public:
    explicit SequenceRecordIterator(std::vector<std::string> &reads_paths);
    ~SequenceRecordIterator();

    bool reset();
    std::optional<GenomeReadData> get_next_record();
    std::vector<ReadFileMetaData> meta;
};


#endif //SRC_SEQUENCERECORDITERATOR_H
