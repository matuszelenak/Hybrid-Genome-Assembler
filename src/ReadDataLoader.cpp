#include "ReadDataLoader.h"

ReadDataLoader::ReadDataLoader(const std::string &path) {
    std::string ext5 = path.substr(path.size() - 5, path.size());
    std::string ext2 = path.substr(path.size() - 2, path.size());
    if (!(ext5 == "fasta" || ext5 == "fastq" || ext2 == "fa" || ext2 == "fq")) {
        throw std::invalid_argument("Only accepts FASTA or FASTQ files");
    }
    if (ext5 == "fasta" || ext2 == "fa") {
        this->read_record = &ReadDataLoader::read_fasta_record;
    } else {
        this->read_record = &ReadDataLoader::read_fastq_record;
    }
    input_file.open(path);
    if (!input_file) {
        throw std::invalid_argument("File does not exist");
    }
}

ReadDataLoader::~ReadDataLoader() {
    if (input_file.is_open()) input_file.close();
}


std::string ReadDataLoader::get_header() {
    if (!input_file.is_open())
        return "";

    std::string header;
    if (!std::getline(input_file, header)) {
        input_file.close();
        return "";
    }
    return header;
}


std::optional<GenomeReadData> ReadDataLoader::read_fasta_record() {
    std::string header = get_header();
    if (header.empty()) return std::nullopt;

    std::string sequence;

    std::getline(input_file, sequence);

    std::vector<Quality> qualities(sequence.length(), 41);

    std::string category = header.substr(0, header.find('-'));
    categories.add_category(category);

    return std::optional<GenomeReadData>{
            {
                    header.substr(1, header.length()),
                    sequence,
                    qualities,
                    category,
                    categories.id_for_category(category)
            }
    };
}

std::optional<GenomeReadData> ReadDataLoader::read_fastq_record() {
    std::string header = get_header();
    if (header.empty()) return std::nullopt;

    std::string sequence, comment, quality;

    std::getline(input_file, sequence);
    std::getline(input_file, comment);
    std::getline(input_file, quality);

    std::vector<Quality> qualities(sequence.length(), -33);

    for (int i = 0; i < sequence.length(); i++) {
        qualities[i] += sequence[i];
    }

    std::string category = header.substr(0, header.find('-'));
    categories.add_category(category);

    return std::optional<GenomeReadData>{
            {
                    header.substr(1, header.length()),
                    sequence,
                    qualities,
                    category,
                    categories.id_for_category(category)
            }
    };
}

std::optional<GenomeReadData> ReadDataLoader::get_next_record() {
    std::lock_guard<std::mutex> lock(_read_mutex);
    return (this->*read_record)();
}