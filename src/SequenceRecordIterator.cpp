#include <fmt/format.h>
#include <iostream>
#include "SequenceRecordIterator.h"
#include "Utils.h"

SequenceRecordIterator::SequenceRecordIterator(std::vector<std::string> &reads_paths) {
    this->paths = reads_paths;
    for (auto &path : reads_paths) {
        total_read_bases.push_back(0);
    }

    current_file_index = 0;
    load_file_at_position(current_file_index);
}

bool SequenceRecordIterator::load_file_at_position(int pos) {
    if (current_file.is_open()){
        current_file.close();
    }

    if (pos >= paths.size()) return false;

    current_file.open(paths[current_file_index]);
    if(!current_file)
    {
        throw std::invalid_argument(fmt::format("File with path \"{}\" does not exist", paths[current_file_index]));
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
    return true;
}

std::string SequenceRecordIterator::get_next_line() {
    std::string in;
    if (!std::getline(current_file, in)) {
        if (load_file_at_position(++current_file_index)){
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

    std::vector<Quality> qualities(sequence.length(), -33);

    for (int i = 0; i < sequence.length(); i++) {
        qualities[i] += sequence[i];
    }

    return {header, sequence, qualities};
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

    show_progress(current_file.tellg(), current_file_size, paths[current_file_index]);

    total_read_bases[current_file_index] += data.sequence.length();

    return std::optional<GenomeReadData>{
            {
                    data.header,
                    data.sequence,
                    data.qualities,
                    (bool) current_file_index
            }
    };
}

SequenceRecordIterator::~SequenceRecordIterator() {
    if (current_file.is_open()) current_file.close();
}
