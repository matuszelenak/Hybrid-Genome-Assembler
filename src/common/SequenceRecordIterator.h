#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <optional>
#include <thread>
#include <mutex>
#include <fmt/format.h>

#include "Types.h"

#ifndef SRC_SEQUENCERECORDITERATOR_H
#define SRC_SEQUENCERECORDITERATOR_H


enum FileType {FASTA, FASTQ};


struct GenomeReadData {
    std::string header = "";
    std::string sequence = "";
    std::string qualities = "";
    CategoryID category_id = 0;

    std::string fastq_string(){
        return fmt::format("{}\n{}\n+\n{}", this->header, this->sequence, this->qualities);
    }
};


struct MetaData {
    uint64_t records = 0;
    uint64_t min_read_length = UINT64_MAX;
    uint64_t max_read_length = 0;
    uint64_t avg_read_length = 0;
    uint64_t total_bases = 0;
};

struct ReadFileMetaData : MetaData {
    std::string filename = "";
    FileType file_type = FASTA;
};


class SequenceRecordIterator {
private:
    std::vector<std::string> paths;
    std::mutex _read_mutex;

    bool _annotate;

    std::ifstream current_file;
    int current_file_index = 0;
    FileType current_file_type = FASTQ;

    uint32_t current_read_index = 1;
    uint32_t show_progress_step = 0;


    std::string get_next_line();

    GenomeReadData read_fastq_record();

    GenomeReadData read_fasta_record();

    GenomeReadData (SequenceRecordIterator::*current_record_method)() = &SequenceRecordIterator::read_fastq_record;

    bool load_file_at_position(int pos);

    void load_meta_data();
public:
    explicit SequenceRecordIterator(std::vector<std::string> &reads_paths, bool annotate);
    ~SequenceRecordIterator();

    bool reset();
    std::optional<GenomeReadData> get_next_record();
    std::vector<ReadFileMetaData> file_meta;
    ReadFileMetaData meta = {};
    uint8_t categories = 1;
};


#endif //SRC_SEQUENCERECORDITERATOR_H
