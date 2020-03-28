#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <optional>
#include <thread>
#include <mutex>

#include "Types.h"

#ifndef SRC_SEQUENCERECORDITERATOR_H
#define SRC_SEQUENCERECORDITERATOR_H


enum FileType {FASTA, FASTQ};


struct GenomeReadData {
    std::string header;
    std::string sequence;
    std::string qualities;
    CategoryID category_id;
};

struct ReadFileMetaData {
    std::string filename;
    FileType file_type;
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
    FileType current_file_type = FASTQ;
    uint64_t current_file_size = 0;


    std::string get_next_line();

    GenomeReadData read_fastq_record();

    GenomeReadData read_fasta_record();

    GenomeReadData (SequenceRecordIterator::*current_record_method)() = &SequenceRecordIterator::read_fastq_record;

    bool load_file_at_position(int pos);

    void get_meta_data();
public:
    explicit SequenceRecordIterator(std::vector<std::string> &reads_paths);
    ~SequenceRecordIterator();

    bool reset();
    std::optional<GenomeReadData> get_next_record();
    std::vector<ReadFileMetaData> meta;

    uint64_t average_read_length();
};


#endif //SRC_SEQUENCERECORDITERATOR_H
