#include <string>
#include <deque>
#include <boost/optional.hpp>
#include <optional>

#include "DNAStructures.h"

#ifndef SRC_KMERITERATOR_H
#define SRC_KMERITERATOR_H


class KmerIterator {
private:
    int kmer_size;
    std::string sequence;

    uint64_t clearing_mask;
    uint8_t complement_shift_by;
    Kmer forward_kmer, complementary_kmer;

    unsigned long position_in_read;

    bool ignore_qualities = false;
    std::vector<Quality > qualities;
    std::deque<Quality> kmer_qualities_window;
    Quality quality_k_behind = 0;
    uint32_t kmer_qualities_sum = 0;

    void roll_forward_strand();
    void roll_complementary_strand();
public:
    explicit KmerIterator(GenomeReadData &read, int k, bool ignore_qualities);

    bool next_kmer();
    std::string number_to_sequence(Kmer kmer);

    Kmer current_kmer = 0;
    KmerQuality current_quality = {0, 0};
};


#endif //SRC_KMERITERATOR_H
