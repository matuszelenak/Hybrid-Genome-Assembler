#include <string>
#include <vector>
#include <fstream>
#include <thread>
#include <mutex>

#ifndef SRC_KMERANALYSISWRITER_H
#define SRC_KMERANALYSISWRITER_H


class KmerAnalysisWriter {
private:
    std::ofstream _out_file;
    std::mutex _read_mutex;
public:
    explicit KmerAnalysisWriter(const std::string &path);
    void write_sequence_data(std::string &sequence_header, std::vector<std::pair<uint64_t, unsigned long> >&positions);
};


#endif //SRC_KMERANALYSISWRITER_H
