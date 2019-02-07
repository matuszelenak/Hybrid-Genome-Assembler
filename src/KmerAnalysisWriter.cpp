#include "KmerAnalysisWriter.h"
#include "KmerIterator.h"

void KmerAnalysisWriter::write_sequence_data(std::string &sequence_header,
                                             std::vector<std::pair<uint64_t, unsigned long> > &positions) {
    std::lock_guard<std::mutex> lock(_read_mutex);
    _out_file << sequence_header << ":[";
    for (auto p: positions){
        _out_file << "(" << p.first << "," << p.second << "),";
    }
    _out_file << "]" << std::endl;
}

KmerAnalysisWriter::KmerAnalysisWriter(const std::string &path) {
    _out_file.open(path, std::ofstream::out | std::ofstream::trunc);
}
