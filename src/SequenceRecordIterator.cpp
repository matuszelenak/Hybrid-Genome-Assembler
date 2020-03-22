#include <fmt/format.h>
#include <iostream>
# include <iomanip>
#include "SequenceRecordIterator.h"
#include "Utils.h"

SequenceRecordIterator::SequenceRecordIterator(std::vector<std::string> &reads_paths) {
    this->paths = reads_paths;
    for (auto &path : reads_paths) {
        auto *file = new std::ifstream(path.c_str());
        if (!file->is_open()) {
            throw std::invalid_argument(fmt::format("File with path \"{}\" does not exist", path));
        }
        this->read_files.push_back(file);
    }
    this->current_record_method = &SequenceRecordIterator::read_fastq_record;
    current_file_index = 0;
    configure_for_file();
}

void SequenceRecordIterator::configure_for_file() {
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

    read_files[current_file_index]->seekg(0, std::ifstream::end);
    current_file_size = read_files[current_file_index]->tellg();
    current_file_position = 0;
    read_files[current_file_index]->seekg(0, std::ifstream::beg);
}

std::string SequenceRecordIterator::get_next_line() {
    std::string in;
    if (!std::getline(*read_files[current_file_index], in)) {
        read_files[current_file_index]->close();
        current_file_index++;
        if (current_file_index < read_files.size()) {
            configure_for_file();
            return get_next_line();
        } else {
            exhausted = true;
            return "";
        }
    }
    current_file_position += (in.length() + 1);
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

    show_progress(current_file_position, current_file_size, paths[current_file_index]);

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
    for (auto file: read_files) {
        if (file->is_open()) file->close();
        free(file);
    }
}
