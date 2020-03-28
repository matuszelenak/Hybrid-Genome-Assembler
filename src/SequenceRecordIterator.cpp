#include <fmt/format.h>
#include <iostream>
#include "SequenceRecordIterator.h"
#include "Utils.h"

SequenceRecordIterator::SequenceRecordIterator(std::vector<std::string> &reads_paths) {
    this->paths = reads_paths;
    get_meta_data();
    reset();
}

bool SequenceRecordIterator::reset() {
    current_file_index = 0;
    bool file_loaded = load_file_at_position(current_file_index);
    if (file_loaded) {
        exhausted = false;
    }
    return file_loaded;
}

void SequenceRecordIterator::get_meta_data() {
    reset();
    meta.resize(paths.size());

    int previous_file_index = -1;
    ReadFileMetaData *current_meta = {};
    std::optional<GenomeReadData> read_record;
    while ((read_record = get_next_record()) != std::nullopt) {
        if (current_file_index != previous_file_index) {
            meta[current_file_index] = {paths[current_file_index].substr(paths[current_file_index].find_last_of("/\\") + 1), 0, UINT64_MAX, 0, 0, 0};
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
        std::cout << fmt::format("{}:\n- {} reads\n- {} total bases\n- {} average read length\n- {} max read length\n- {} min read length\n\n",
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
        if (!comment.empty() && comment[0] == '+')
            this->current_record_method = &SequenceRecordIterator::read_fastq_record;
    } else if (header[0] == '>') {
        this->current_record_method = &SequenceRecordIterator::read_fasta_record;
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
            exhausted = true;
            return "";
        }
    }
    return in;
}

SeqRecordData SequenceRecordIterator::read_fastq_record(std::string &header) {
    std::string sequence, comment, quality;
    sequence = get_next_line();
    comment = get_next_line();
    quality = get_next_line();

//    std::vector<Quality> qualities(sequence.length(), -33);
//
//    for (int i = 0; i < sequence.length(); i++) {
//        qualities[i] += sequence[i];
//    }

    return {header, sequence, {}};
}

SeqRecordData SequenceRecordIterator::read_fasta_record(std::string &header) {
    std::string sequence;
    sequence = get_next_line();

    std::vector<Quality> qualities(sequence.length(), 40);

    return {header, sequence, qualities};
}

std::optional<GenomeReadData> SequenceRecordIterator::get_next_record() {
    std::lock_guard<std::mutex> lock(_read_mutex);
    if (exhausted) return std::nullopt;

    std::string header = get_next_line();
    if (header.empty()) return std::nullopt;

    SeqRecordData data = (this->*current_record_method)(header);

    if (data.header.empty() || data.sequence.empty()) {
        return std::nullopt;
    }

    //show_progress(current_file.tellg(), current_file_size, paths[current_file_index]);

    return std::optional<GenomeReadData>{
            {
                    data.header,
                    data.sequence,
                    data.qualities,
                    current_file_index
            }
    };
}

SequenceRecordIterator::~SequenceRecordIterator() {
    if (current_file.is_open()) current_file.close();
}