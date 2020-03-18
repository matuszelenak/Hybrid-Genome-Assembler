#include "SequenceReader.h"

SequenceReader::SequenceReader(const std::string &path) {
    std::string ext5 = path.substr(path.size() - 5, path.size());
    std::string ext2 = path.substr(path.size() - 2, path.size());
    if (!(ext5 == "fasta" || ext5 == "fastq" || ext2 == "fa" || ext2 == "fq")) {
        throw std::invalid_argument("Only accepts FASTA or FASTQ files");
    }
    if (ext5 == "fasta" || ext2 == "fa") {
        this->read_sequence_line = &SequenceReader::read_fasta_sequence;
    } else {
        this->read_sequence_line = &SequenceReader::read_fastq_sequence;
    }
    input_file.open(path);
    if (!input_file) {
        throw std::invalid_argument("File does not exist");
    }
}

std::optional<GenomeRead> SequenceReader::read_fasta_sequence() {
    if (!input_file.is_open())
        return std::nullopt;

    std::string header, sequence;
    if (!std::getline(input_file, header)) {
        input_file.close();
        return std::nullopt;
    }
    std::getline(input_file, sequence);
    std::vector<Quality> qualities;

    return std::optional<GenomeRead>{{header.substr(1, header.length()), sequence, qualities, header.substr(1, header.find(' ') - 1)}};
}

std::optional<GenomeRead> SequenceReader::read_fastq_sequence() {
    if (!input_file.is_open())
        return std::nullopt;

    std::string header, sequence, comment, quality, category;
    if (!std::getline(input_file, header)) {
        input_file.close();
        return std::nullopt;
    }
    std::getline(input_file, sequence);
    std::getline(input_file, comment);
    std::getline(input_file, quality);

    std::vector<Quality> qualities;
    for (char q: quality) {
        qualities.push_back(q - 33);
    }

    return std::optional<GenomeRead>{{header.substr(1, header.length()), sequence, qualities, header.substr(header.find(' '), header.length() - header.find(' '))}};
}

std::optional<GenomeRead> SequenceReader::get_next_record() {
    std::lock_guard<std::mutex> lock(_read_mutex);
    return (this->*read_sequence_line)();
}
