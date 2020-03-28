#include <fmt/format.h>
#include <iostream>
#include <numeric>

#include "SequenceRecordIterator.h"
#include "Utils.h"

SequenceRecordIterator::SequenceRecordIterator(std::vector<std::string> &reads_paths) {
    this->paths = reads_paths;
    get_meta_data();
    reset();
}

bool SequenceRecordIterator::reset() {
    current_file_index = 0;
    return load_file_at_position(current_file_index);
}

void SequenceRecordIterator::get_meta_data() {
    reset();
    meta.resize(paths.size());

    int previous_file_index = -1;
    ReadFileMetaData *current_meta = {};
    std::optional<GenomeReadData> read_record;
    while ((read_record = get_next_record()) != std::nullopt) {
        if (current_file_index != previous_file_index) {
            meta[current_file_index] = {paths[current_file_index].substr(paths[current_file_index].find_last_of("/\\") + 1), current_file_type, 0, UINT64_MAX, 0, 0, 0};
            current_meta = &meta[current_file_index];
            previous_file_index = current_file_index;
        }

        uint64_t seq_length = read_record->sequence.length();

        current_meta->total_bases += seq_length;
        current_meta->min_read_length = std::min(current_meta->min_read_length, seq_length);
        current_meta->max_read_length = std::max(current_meta->max_read_length, seq_length);
        current_meta->records++;
        current_meta->avg_read_length += seq_length;
    }
    for (auto &data : meta) {
        data.avg_read_length /= data.records;
        std::cout << fmt::format("{}:\n-{} reads\n- {} total bases\n- {} average read length\n- {} max read length\n- {} min read length\n\n",
                                 data.filename, data.records, data.total_bases, data.avg_read_length, data.max_read_length, data.min_read_length);
    }
}

bool SequenceRecordIterator::load_file_at_position(int pos) {
    if (current_file.is_open()) {
        current_file.close();
    }

    if (pos >= paths.size()) return false;

    current_file.open(paths[pos]);
    if (!current_file) {
        throw std::invalid_argument(fmt::format("File with path \"{}\" does not exist", paths[pos]));
    }

    //Measure the size of the file
    current_file.seekg(0, std::ifstream::end);
    current_file_size = current_file.tellg();
    current_file.seekg(0, std::ifstream::beg);

    // Determine the file format
    std::string header = get_next_line();
    std::string sequence = get_next_line();
    if (header.empty()) {
        throw std::logic_error("File is empty");
    }

    if (header[0] == '@') {
        std::string comment = get_next_line();
        if (!comment.empty() && comment[0] == '+') {
            this->current_record_method = &SequenceRecordIterator::read_fastq_record;
            this->current_file_type = FASTQ;
        }
    } else if (header[0] == '>') {
        this->current_record_method = &SequenceRecordIterator::read_fasta_record;
        this->current_file_type = FASTA;
    } else
        throw std::logic_error("Unrecognized file format");

    current_file.seekg(0, std::ifstream::beg);
    current_file_index = pos;
    return true;
}

std::string SequenceRecordIterator::get_next_line() {
    std::string in;
    if (!std::getline(current_file, in)) {
        if (load_file_at_position(current_file_index + 1)) {
            return get_next_line();
        } else {
            throw std::length_error("");
        }
    }
    return in;
}

GenomeReadData SequenceRecordIterator::read_fastq_record() {
    std::string header, sequence, comment, qualities;
    header = get_next_line();
    sequence = get_next_line();
    comment = get_next_line();
    qualities = get_next_line();

    return {header, sequence, qualities, current_file_index};
}

GenomeReadData SequenceRecordIterator::read_fasta_record() {
    std::string header, sequence;
    header = get_next_line();
    sequence = get_next_line();

    return {header, sequence, "", current_file_index};
}

std::optional<GenomeReadData> SequenceRecordIterator::get_next_record() {
    std::lock_guard<std::mutex> lock(_read_mutex);
    try {
        return std::optional<GenomeReadData>{(this->*current_record_method)()};
    } catch (const std::length_error &e){
        return std::nullopt;
    }
}

SequenceRecordIterator::~SequenceRecordIterator() {
    if (current_file.is_open()) current_file.close();
}

uint64_t SequenceRecordIterator::average_read_length() {
    return std::accumulate(meta.begin(), meta.end(), 0, [](uint64_t acc, ReadFileMetaData &m) -> uint64_t { return acc + m.avg_read_length; }) / meta.size();
}
