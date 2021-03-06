#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <optional>
#include <thread>
#include <mutex>
#include <fmt/format.h>
#include <boost/regex.hpp>

#include "Types.h"

#ifndef SRC_SEQUENCERECORDITERATOR_H
#define SRC_SEQUENCERECORDITERATOR_H

typedef uint32_t ReadID;


enum FileType {
    FASTA, FASTQ, UNKNOWN
};


struct GenomeReadData {
    ReadID id = 0;
    std::string header = "";
    std::string sequence = "";
    std::string qualities = "";

    //DEBUG
    CategoryID category_id = 0;
    uint32_t start = 0;
    uint32_t end = 0;
    //END DEBUG

    std::string fasta_string(){
        return fmt::format(">{}\n{}", this->header, this->sequence);
    }

    std::string fastq_string() {
        return fmt::format("@{}\n{}\n+\n{}", this->header, this->sequence, this->qualities);
    }

    std::string fastX_string() {
        if (this->qualities.empty()){
            return this->fasta_string();
        } else return this->fastq_string();
    }
};


struct MetaData {
    std::string filename = "";
    uint64_t records = 0;
    uint64_t min_read_length = UINT64_MAX;
    uint64_t max_read_length = 0;
    uint64_t avg_read_length = 0;
    uint64_t total_bases = 0;

    std::string repr() {
        return fmt::format("{}:\n- {} reads\n- {} total bases\n- {} average read length\n- {} max read length\n- {} min read length\n\n",
                           this->filename, this->records, this->total_bases, this->avg_read_length, this->max_read_length, this->min_read_length);
    }
};

struct ReadFileMetaData : MetaData {
    FileType file_type = UNKNOWN;
};


class SequenceRecordIterator {
private:
    std::vector<std::string> paths;
    std::mutex _read_mutex;

    bool _annotate = true;

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

    boost::regex simlord_header_regex{";length=([0-9]+)bp;startpos=([0-9]+);"};
    std::pair<uint32_t, uint32_t> parse_simlord_header(std::string &header);

    boost::regex nanosimh_header_regex{"_([0-9]+)_[^_]+_[^_]+_[^_]+_[^_]+_([0-9]+)_"};
    std::pair<uint32_t, uint32_t> parse_nanosimh_header(std::string &header);

    boost::regex PaSS_header_regex{"([0-9]+)_([0-9]+)\\|([0-9]+)\\|"};
    std::pair<uint32_t, uint32_t> parse_PaSS_header(std::string &header);

    std::pair<uint32_t, uint32_t> parse_unknown_header(std::string &header){ return {0, 0}; };

    std::pair<uint32_t, uint32_t> (SequenceRecordIterator::*parse_header_method)(std::string &header) = &SequenceRecordIterator::parse_unknown_header;

public:
    explicit SequenceRecordIterator(std::string &path);

    SequenceRecordIterator(std::vector<std::string> &reads_paths, bool annotate);

    ~SequenceRecordIterator();

    bool rewind();

    std::optional<GenomeReadData> get_next_record();

    std::vector<ReadFileMetaData> file_meta;
    ReadFileMetaData meta;
    uint8_t categories = 1;

    bool show_progress = false;
};


#endif //SRC_SEQUENCERECORDITERATOR_H
