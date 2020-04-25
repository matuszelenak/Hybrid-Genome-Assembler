#include <fmt/format.h>
#include <iostream>
#include <boost/algorithm/string/join.hpp>

#include "SequenceRecordIterator.h"
#include "Utils.h"


SequenceRecordIterator::SequenceRecordIterator(std::string &path) {
    paths = {path};
    load_meta_data();
    rewind();
    show_progress = true;
}


SequenceRecordIterator::SequenceRecordIterator(std::vector<std::string> &reads_paths, bool annotate) {
    _annotate = annotate;
    this->paths = reads_paths;
    load_meta_data();
    categories = annotate ? file_meta.size() : 1;
    rewind();
    show_progress = true;
}

bool SequenceRecordIterator::rewind() {
    current_file_index = 0;
    current_read_index = 1;
    return load_file_at_position(current_file_index);
}

void SequenceRecordIterator::load_meta_data() {
    rewind();
    file_meta.resize(paths.size());

    int previous_file_index = -1;
    std::vector<std::string> filenames;
    ReadFileMetaData *current_meta = {};
    std::optional<GenomeReadData> read_record;
    while ((read_record = get_next_record()) != std::nullopt) {
        if (current_file_index != previous_file_index) {
            file_meta[current_file_index] = {};
            file_meta[current_file_index].filename = paths[current_file_index].substr(paths[current_file_index].find_last_of("/\\") + 1);
            file_meta[current_file_index].file_type = current_file_type;

            current_meta = &file_meta[current_file_index];
            previous_file_index = current_file_index;

            filenames.push_back(current_meta->filename);
        }

        uint64_t seq_length = read_record->sequence.length();

        current_meta->total_bases += seq_length;
        current_meta->min_read_length = std::min(current_meta->min_read_length, seq_length);
        current_meta->max_read_length = std::max(current_meta->max_read_length, seq_length);
        current_meta->records++;
        current_meta->avg_read_length += seq_length;

        meta.total_bases += seq_length;
        meta.min_read_length = std::min(meta.min_read_length, seq_length);
        meta.records++;
        meta.avg_read_length += seq_length;
    }
    meta.avg_read_length /= meta.records;
    meta.filename = boost::algorithm::join(filenames, "__");
    show_progress_step = std::max(1ul, meta.records / 1000);

    for (auto &data : file_meta) {
        data.avg_read_length /= data.records;
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

    // Determine the file format
    try{
        std::string header = get_next_line();
        std::string sequence = get_next_line();

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

        auto position_and_length = parse_simlord_header(header);
        if (position_and_length.second != 0){
            this->parse_header_method = &SequenceRecordIterator::parse_simlord_header;
        }

        position_and_length = parse_nanosimh_header(header);
        if (position_and_length.second != 0){
            this->parse_header_method = &SequenceRecordIterator::parse_nanosimh_header;
        }

        position_and_length = parse_PaSS_header(header);
        if (position_and_length.second != 0){
            this->parse_header_method = &SequenceRecordIterator::parse_PaSS_header;
        }

    } catch (const std::length_error &e) {
        throw std::logic_error("File is empty");
    }

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

    return {current_read_index, header, sequence, qualities, _annotate ? current_file_index : 0, 0, 0};
}

GenomeReadData SequenceRecordIterator::read_fasta_record() {
    std::string header, sequence;
    header = get_next_line();
    sequence = get_next_line();

    return {current_read_index, header, sequence, "", _annotate ? current_file_index : 0, 0, 0};
}

std::optional<GenomeReadData> SequenceRecordIterator::get_next_record() {
    std::lock_guard<std::mutex> lock(_read_mutex);
    try {
        auto result = std::optional<GenomeReadData>{(this->*current_record_method)()};
        auto position_and_length = (this->*parse_header_method)(result->header);
        if (position_and_length.first != 0){
            result->start = position_and_length.first;
            result->end = position_and_length.first + position_and_length.second;
        }

        if (show_progress && show_progress_step && (current_read_index % show_progress_step == 0 || current_read_index == meta.records)){
            progress_bar(current_read_index, meta.records, "Iterating reads...");
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

std::pair<uint32_t, uint32_t> SequenceRecordIterator::parse_simlord_header(std::string &header) {
    boost::smatch match;
    if (boost::regex_search(header, match, simlord_header_regex))
    {
        return {std::stoul(match[2].str()), std::stoul(match[1].str())};
    }
    return {0, 0};
}

std::pair<uint32_t, uint32_t> SequenceRecordIterator::parse_nanosimh_header(std::string &header) {
    boost::smatch match;
    if (boost::regex_search(header, match, nanosimh_header_regex))
    {
        return {std::stoul(match[1].str()), std::stoul(match[2].str())};
    }
    return {0, 0};
}

std::pair<uint32_t, uint32_t> SequenceRecordIterator::parse_PaSS_header(std::string &header) {
    boost::smatch match;
    if (boost::regex_search(header, match, PaSS_header_regex))
    {
        auto length = std::stoul(match[2].str()) - std::stoul(match[1].str());
        return {std::stoul(match[3].str()), length};
    }
    return {0, 0};
}
