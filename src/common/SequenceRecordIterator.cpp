#include <fmt/format.h>
#include <iostream>
#include <numeric>

#include "SequenceRecordIterator.h"
#include "Utils.h"

SequenceRecordIterator::SequenceRecordIterator(std::vector<std::string> &reads_paths) {
    this->paths = reads_paths;
    load_meta_data();
    reset();
}

bool SequenceRecordIterator::reset() {
    current_file_index = 0;
    current_read_index = 1;
    return load_file_at_position(current_file_index);
}

void SequenceRecordIterator::load_meta_data() {
    reset();
    file_meta.resize(paths.size());

    int previous_file_index = -1;
    ReadFileMetaData *current_meta = {};
    std::optional<GenomeReadData> read_record;
    while ((read_record = get_next_record()) != std::nullopt) {
        if (current_file_index != previous_file_index) {
            file_meta[current_file_index] = {};
            file_meta[current_file_index].filename = paths[current_file_index].substr(paths[current_file_index].find_last_of("/\\") + 1);
            file_meta[current_file_index].file_type = current_file_type;

            current_meta = &file_meta[current_file_index];
            previous_file_index = current_file_index;
        }

        uint64_t seq_length = read_record->sequence.length();

        current_meta->total_bases += seq_length;
        current_meta->min_read_length = std::min(current_meta->min_read_length, seq_length);
        current_meta->max_read_length = std::max(current_meta->max_read_length, seq_length);
        current_meta->records++;
        current_meta->avg_read_length += seq_length;

        meta.total_bases += seq_length;
        meta.min_read_length = std::min(meta.min_read_length, seq_length);
        meta.max_read_length = std::max(meta.max_read_length, seq_length);
        meta.records++;
        meta.avg_read_length += seq_length;
    }
    meta.avg_read_length /= meta.records;
    show_progress_step = std::max(1ul, meta.records / 1000);

    for (auto &data : file_meta) {
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
        auto result = std::optional<GenomeReadData>{(this->*current_record_method)()};
        if (show_progress_step && (current_read_index % show_progress_step == 0 || current_read_index == meta.records)){
            show_progress(current_read_index, meta.records, "Iterating reads...");
        }
        current_read_index++;
        return result;
    } catch (const std::length_error &e){
        return std::nullopt;
    }
}

SequenceRecordIterator::~SequenceRecordIterator() {
    if (current_file.is_open()) current_file.close();
}

uint64_t SequenceRecordIterator::average_read_length() {
    return std::accumulate(file_meta.begin(), file_meta.end(), 0, [](uint64_t acc, ReadFileMetaData &m) -> uint64_t { return acc + m.avg_read_length; }) / file_meta.size();
}
