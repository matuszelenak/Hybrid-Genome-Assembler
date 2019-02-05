//
// Created by whiskas on 04/02/19.
//

#include "SequenceReader.h"

SequenceReader::SequenceReader(const std::string &path){
    std::string ext5 = path.substr(path.size() - 5, path.size());
    std::string ext2 = path.substr(path.size() - 2, path.size());
    if (!(ext5 == "fasta" || ext5 == "fastq" || ext2 == "fa" || ext2 == "fq")){
        throw std::invalid_argument("Only accepts FASTA or FASTQ files");
    }
    if (ext5 == "fasta" || ext2 == "fa"){
        this->read_sequence_line = &SequenceReader::read_fasta_sequence;
    }
    else{
        this->read_sequence_line = &SequenceReader::read_fastq_sequence;
    }
    input_file.open(path);
    if(!input_file)
    {
        throw std::invalid_argument("File does not exist");
    }
}

std::optional<std::string> SequenceReader::read_fasta_sequence() {
    if (!input_file.is_open()) return "";
    std::string header, sequence;
    if (!(input_file >> header)) {
        input_file.close();
        return std::nullopt;
    }
    input_file >> sequence;
    return sequence.substr(0, sequence.size() - 1);
}

std::optional<std::string> SequenceReader::read_fastq_sequence() {
    if (!input_file.is_open()) return "";
    std::string line, sequence;
    if (!(input_file >> line)){
        input_file.close();
        return std::nullopt;
    }
    input_file >> sequence;
    input_file >> line; // + and comment
    input_file >> line; // sequence quality, not used for now
    return sequence.substr(0, sequence.size() - 1);
}

std::optional<std::string> SequenceReader::get_next_record(){
    return (this->*read_sequence_line)();
}